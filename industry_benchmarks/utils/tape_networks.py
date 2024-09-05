import json
import click
import pathlib
from glob import glob
import networkx as nx

from gufe.tokenization import JSON_HANDLER
import openfe
from openfe import ChemicalSystem, LigandAtomMapping, Protocol
from openfe import Transformation, AlchemicalNetwork, LigandNetwork
from openfe import SolventComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe.setup import KartografAtomMapper, lomap_scorers
from konnektor.network_planners import MstConcatenator
from konnektor.network_tools import concatenate_networks
import plan_rbfe_network

try:
    from konnektor.network_planners import MstConcatenator
    from konnektor.network_tools import concatenate_networks
except ModuleNotFoundError:  # This will attempt to install kartograf
    print("Did not find Konnektor! will attempt install!")
    os.system(
        "pip install git+https://github.com/OpenFreeEnergy/konnektor.git"
    )
    from konnektor.network_planners import MstConcatenator
    from konnektor.network_tools import concatenate_networks


def parse_alchemical_network(
    input_alchem_network_json_path: str,
) -> AlchemicalNetwork:
    j_dict = json.load(
        open(input_alchem_network_json_path, "r"), cls=JSON_HANDLER.decoder
    )
    alchem_network = AlchemicalNetwork.from_dict(j_dict)
    return alchem_network


def parse_result_folders(
    result_files_regex: str, input_ligand_network: LigandNetwork = None
) -> AlchemicalNetwork:
    nodes = []
    edges = {}
    edge_count = {}
    for resp in glob(result_files_regex):
        res =  json.load(open(resp, 'rb'), cls=JSON_HANDLER.decoder)
        if any("exception" in u for u in res["unit_results"].values()):
            continue

        for k in res["protocol_result"]["data"]:
            units = res["protocol_result"]["data"][k]

            name = units[0]["name"].split(" repeat")[0]
            if "stateA" in units[0]["inputs"]:
                stateA = ChemicalSystem.from_dict(
                    units[0]["inputs"]["stateA"]
                )
                stateB = ChemicalSystem.from_dict(
                    units[0]["inputs"]["stateB"]
                )
                ligmap = LigandAtomMapping.from_dict(
                    units[0]["inputs"]["ligandmapping"]
                )
            else:  # inputs were cleaned out! reengineer - hacky alchem net.
                if input_ligand_network is not None:

                    transfromation_name = (
                        units[0]["name"].split(" repeat")[0].strip()
                    )
                    componentA_name, componentB_name = (
                        transfromation_name.split(" to ")
                    )

                    ligmap = None
                    for e in input_ligand_network.edges:
                        if all(
                            cn in (e.componentA.name, e.componentB.name)
                            for cn in [componentA_name, componentB_name]
                        ):
                            ligmap = e
                            break
                    if ligmap is None:

                        raise ValueError(
                            "could not finde edge: ("
                            + e.componentA.name
                            + ", "
                            + e.componentB.name
                            + ")"
                        )
                    else:
                        if "solvent" in str(units[0]["outputs"]["nc"]):
                            stateA = ChemicalSystem(
                                components={
                                    "ligand": ligmap.componentA,
                                    "solvent": SolventComponent(),
                                }
                            )
                            stateB = ChemicalSystem(
                                components={
                                    "ligand": ligmap.componentB,
                                    "solvent": SolventComponent(),
                                }
                            )
                        else:
                            stateA = ChemicalSystem(
                                components={"ligand": ligmap.componentA}
                            )
                            stateB = ChemicalSystem(
                                components={"ligand": ligmap.componentB}
                            )
                else:
                    raise ValueError(
                        "Can not read alchemical input, as components are missing, please provied input ligand_network."
                    )

            # Count if both transformations present.
            if name not in edge_count:
              edge_count[name] = {"protein": 0, "solvent":0}
                
            if (
                "protein" in stateA.components
                or "solvent" not in stateA.components
              ):
              edge_count[name]["protein"] += 1
            else:
              edge_count[name]["solvent"] += 1
          
            try:
                protocol = Protocol.from_dict(units[0]["inputs"]["protocol"])
            except:
                protocol = None
            t = Transformation(
                stateA=stateA,
                stateB=stateB,
                mapping=ligmap,
                name=name,
                protocol=protocol,
            )

            nodes.extend([stateA, stateB])
            if name not in edges:
                edges[name] = t

    # get only finished edges:
    fin_edge = []
    for name in edges:
        edge = edges[name]
        ec = edge_count[name]

        if ec["solvent"] >= 1 and ec["protein"] >= 1:
            fin_edge.append(edge)

    alchemical_network = AlchemicalNetwork(nodes=nodes, edges=fin_edge)
    return alchemical_network


def alchemical_network_to_ligand_network(alchemical_network) -> LigandNetwork:
    edges = []
    nodes = []
    for e in alchemical_network.edges:
        nA = e.stateA.components["ligand"]
        nB = e.stateB.components["ligand"]

        edges.append(e.mapping)
        nodes.append(nA)
        nodes.append(nB)

    network = LigandNetwork(nodes=set(nodes), edges=set(edges))
    return network


def ligand_network_to_networkx(ligand_network) -> nx.Graph:

    edges = ligand_network.edges
    wedges = []
    for edge in list(edges):
        wedges.append(
            [edge.componentA, edge.componentB, edge.annotations["score"]]
        )

    g = nx.Graph(nodes=ligand_network.nodes)
    g.add_weighted_edges_from(ebunch_to_add=wedges)
    return g



def decomposite_disconnected_ligand_network(
    ligand_network: LigandNetwork,
) -> list[LigandNetwork]:
    # build alchemical sub-networks
    g = ligand_network_to_networkx(ligand_network)
    all_edges = list(ligand_network.edges)
    sub_network_nodes = list(nx.connected_components(g))

    ligand_sub_networks = []
    for i, snn in enumerate(sub_network_nodes):
        # collect edges
        sub_edges = []
        for e in all_edges:
            if any(n in (e.componentA, e.componentB) for n in snn):
                sub_edges.append(e)

        ligand_sub_network = LigandNetwork(nodes=snn, edges=sub_edges)
        ligand_sub_networks.append(ligand_sub_network)

    return ligand_sub_networks


def get_new_network_tapes(
    ligand_sub_networks: list[LigandNetwork], 
    input_ligand_network: LigandNetwork,
    n_connecting_edges: int = 3

) -> LigandNetwork:

    mapper = KartografAtomMapper()
    scorer = lomap_scorers.default_lomap_score
    concatenator = MstConcatenator(
        mapper=mapper,
        scorer=scorer,
        n_connecting_edges=n_connecting_edges,
        n_processes=1,
    )

    in_edges = input_ligand_network.edges

    concatenated_network = ligand_sub_networks[0]
    tape_edges = []
    for ligand_sub_network in ligand_sub_networks[1:]:
        nedges= min(n_connecting_edges*10, len(concatenated_network.edges))
        
        concatenator.n_connecting_edges = nedges
        concatenated_network = concatenate_networks(
            [ligand_sub_network, concatenated_network],
            concatenator=concatenator,
        )
        concatenated_edges = concatenated_network.edges
        possible_edges = list(sorted(list(concatenated_edges.difference(in_edges)), key=lambda e: e.annotations["score"], reverse=True))

        if len(possible_edges) >= n_connecting_edges:
            tape_edges.extend(possible_edges[:n_connecting_edges])
        else:
             tape_edges.extend(possible_edges)
    tape_nodes = set(
        [n for e in tape_edges for n in [e.componentA, e.componentB]]
    )
    return LigandNetwork(nodes=tape_nodes, edges=tape_edges)


def get_taped_alchemical_network(ducktape_network, alchemical_network):
    solv = openfe.SolventComponent()

    for node in alchemical_network.nodes:
        if "protein" in node.components:
            prot = node.components["protein"]
            # Add cofactors if present
            cofactors = []
            if len(node.components) > 3:
                number_cofactors = len(node.components) - 3
                for i in range(number_cofactors):
                    cofactor_name = f"cofactor_{i}"
                    cofactors.append(node.components[cofactor_name])
            break

    # Create the AlchemicalTransformations, and storing them to an
    # AlchemicalNetwork
    transformations = []
    for mapping in ducktape_network.edges:
        # Get different settings depending on whether the transformation
        # involves a change in net charge
        charge_difference = plan_rbfe_network.get_alchemical_charge_difference(mapping)
        if abs(charge_difference) > 1e-3:
            # Raise a warning that a charge changing transformation is included
            # in the network
            wmsg = ("Charge changing transformation between ligands "
                    f"{mapping.componentA.name} and "
                    f"{mapping.componentB.name}. "
                    "A more expensive protocol with 22 lambda windows, "
                    "sampled "
                    "for 20 ns each, will be used here.")
            warnings.warn(wmsg)
            # Get settings for charge changing transformations
            rfe_settings = plan_rbfe_network.get_settings_charge_changes()
        else:
            rfe_settings = plan_rbfe_network.get_settings()
        for leg in ['solvent', 'complex']:
            # use the solvent and protein created above
            sysA_dict = {'ligand': mapping.componentA,
                         'solvent': solv}
            sysB_dict = {'ligand': mapping.componentB,
                         'solvent': solv}

            if leg == 'complex':
                sysA_dict['protein'] = prot
                sysB_dict['protein'] = prot
                if len(cofactors) > 0:

                    for cofactor, entry in zip(cofactors_smc,
                                               string.ascii_lowercase):
                        cofactor_name = f"cofactor_{entry}"
                        sysA_dict[cofactor_name] = cofactor
                        sysB_dict[cofactor_name] = cofactor

            sysA = openfe.ChemicalSystem(sysA_dict)
            sysB = openfe.ChemicalSystem(sysB_dict)

            name = (f"{leg}_{mapping.componentA.name}_"
                    f"{mapping.componentB.name}")

            rbfe_protocol = RelativeHybridTopologyProtocol(
                settings=rfe_settings)
            transformation = openfe.Transformation(
                stateA=sysA,
                stateB=sysB,
                mapping=mapping,
                protocol=rbfe_protocol,
                name=name,
            )

            transformations.append(transformation)

    # Create the taped alchemical network
    alchemical_network = openfe.AlchemicalNetwork(transformations)

    return alchemical_network


def do_taping(
    result_files_regex: str,
    input_alchem_network_file: str,
    output_alchemical_network_folder: str = "./ducktaping_transformations",
    n_connecting_edges: int = 3,
):
    print("#"*40)
    print("Network Taping")
    print("#"*40)
    print()
    print("Parsing input files:")
    # Parse Ligand Network
    print(f"\treading input Alchemical Network")
    input_alchem_network = parse_alchemical_network(input_alchem_network_file)
    input_ligand_network = alchemical_network_to_ligand_network(input_alchem_network)
    
    # Parse Alchemical Network
    print(f"\treading result files")
    res_alchemical_network = parse_result_folders(result_files_regex, input_ligand_network)
    res_ligand_network = alchemical_network_to_ligand_network(res_alchemical_network)

    
    # decomposite disconnected sets
    ligand_sub_networks = decomposite_disconnected_ligand_network(res_ligand_network)
    
    # Add missing nodes to ligand_sub_networks if ligand network provided
    # get network node differences
    alchemical_network_ligands = set([n for n in res_ligand_network.nodes])
    missing_commponents = list(input_ligand_network.nodes.difference(res_ligand_network.nodes))
    
    for missing_commponent in missing_commponents:
        ligand_sub_networks.append(LigandNetwork(nodes=[missing_commponent], edges=[]))
    ligand_sub_networks = list(sorted(ligand_sub_networks, key=lambda sn: len(sn.nodes), reverse=True))
    
    # check status
    print(f"Inspecting input files:")
    print(f"\tinput LigandNetwork  nodes: {len(input_ligand_network.nodes)}")
    print(f"\tinput LigandNetwork  edges: {len(input_ligand_network.edges)}")    
    print(f"\tresult LigandNetwork nodes: {len(res_ligand_network.nodes)}")
    print(f"\tresult LigandNetwork edges: {len(res_ligand_network.edges)}")
    print()
    
    print(f"Found the following problems:")
    print(f"\tmissing components: {len(missing_commponents)}")
    print(f"\tdisconnected networks: {len(ligand_sub_networks)}")
    print(f"\tdisconnected networks sizes: {[len(sn.nodes) for sn in ligand_sub_networks]}")
    
    if len(ligand_sub_networks) == 1:
        raise ValueError("Did not find disconnected components in alchemical network!")

    print()
    
    # Tape the networks together
    print("Taping the Networks together:")
    network_tapes = get_new_network_tapes(ligand_sub_networks=ligand_sub_networks, input_ligand_network=input_ligand_network, n_connecting_edges = 3)
    
    print(f"\tGenerated new taping LigandMapping edges: {len(network_tapes.edges)}")
    print(f"\tUsing ligand nodes for taping: {len(network_tapes.nodes)}")
    print()
    
    # write out
    print("Write out the tapes for the network:")
    output_alchemical_network_folder.mkdir(exist_ok=False, parents=True)
    duck_taped_alchemical_network = get_taped_alchemical_network(duck_taped_network, input_alchem_network)
    alchemical_network_json_fp = output_alchemical_network_folder / "alchemical_network.json"
    json.dump(
        duck_taped_alchemical_network.to_dict(),
        alchemical_network_json_fp.open(mode="w"),
        cls=JSON_HANDLER.encoder
    )

    # Write out each transformation
    # Create a subdirectory for the transformations
    transforms_dir = pathlib.Path(output_alchemical_network_folder / "transformations")
    transforms_dir.mkdir(exist_ok=True, parents=True)

    for transform in duck_taped_alchemical_network.edges:
        transform.dump(transforms_dir / f"{transform.name}.json")

    print("Done!")


@click.command
@click.option(
    "-rp",
    "--results_regex_path",
    type=click.Path(dir_okay=False, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("./results*/*json"),
    help="regular expression to find all result files of this system. (default: 'results*/*json')",
)
@click.option(
    "-ia",
    "--input_alchem_network_file",
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./alchemicalNetwork/alchemical_network.json"),
    help="Path to the used ligand network plan (default: 'alchemicalNetwork/ligand_network.graphml')",
)
@click.option(
    "-oa",
    "--output_alchemical_network",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("./tapesForAlchemicalNetwork"),
    help="Directory name in which to store the additional taping transformation json files. (default: './tapesForAlchemicalNetwork')",
)

@click.option(
    "-ne",
    "--n_connecting_edges",
    type=int,
    default=3,
    help="number of edges per connection being established.",
)
def cli_do_taping(
    results_regex_path: str,
    input_alchem_network_file: str,
    output_alchemical_network: str,
    n_connecting_edges: int,
):
    do_taping(
        result_files_regex=str(results_regex_path),
        input_alchem_network_file=input_alchem_network_file,
        output_alchemical_network_folder=output_alchemical_network,
        n_connecting_edges=n_connecting_edges,
    )


# main Runs
if __name__ == "__main__":
    cli_do_taping()
