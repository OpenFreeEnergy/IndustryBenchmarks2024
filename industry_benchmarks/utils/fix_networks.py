import json
import os
import pathlib
import networkx as nx
import warnings

import gufe
from gufe.tokenization import JSON_HANDLER
import openfe
from openfe import ChemicalSystem, LigandAtomMapping
from openfe import Transformation, AlchemicalNetwork, LigandNetwork
from openfe import ProteinComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import (
    RelativeHybridTopologyProtocol
)

from openfe.setup import KartografAtomMapper, lomap_scorers


def parse_alchemical_network(
    input_alchem_network_json_path: str,
) -> AlchemicalNetwork:
    """
    Read an AlchemicalNetwork
    """
    j_dict = json.load(
        open(input_alchem_network_json_path, "r"), cls=JSON_HANDLER.decoder
    )
    alchem_network = AlchemicalNetwork.from_dict(j_dict)
    return alchem_network


def _get_check_results_json(filename: str):
    """
    Get a json file, open it, check it's
    a results json, and return it.
    """
    if not os.path.exists(filename):
        print(f"Error: {filename} does not exist.")
        return None

    ru = json.load(open(filename, 'rb'), cls=JSON_HANDLER.decoder)

    if ru.get("__qualname__") == "AlchemicalNetwork":
        print(f"{filename} is a network_setup.json, skipping")
        return None
    if ru.get("__qualname__") == "Transformation":
        print(f"{filename} is an input json, skipping")
        return None
    if "unit_results" not in ru.keys():
        print(f"{filename} has no unit results - likely a failure")
        return None
    if any('exception' in u for u in ru['unit_results'].values()):
        print(f"{filename} is a failed simulation")
        return None

    return ru


def get_transformation_alternate(
    pur: gufe.ProtocolUnitResult, ligand_network: LigandNetwork
):
    """
    Getting a transformation if things got deleted.
    """
    # The way this was being done in PR #144 doesn't work because
    # it would lead to either a solvent or vacuum simulation, not a
    # complex one!
    errmsg = ("The transformation data is not available, "
              "please contact the OpenFE team")
    raise ValueError(errmsg)


def get_transformation(pur: gufe.ProtocolUnitResult):
    """
    Get a transformation out of a PUR

    Returns
    -------
    transform : gufe.Transformation
      The transformation for this edge
    phase : str
      Either "complex" or "solvent" depending on the phase type.
    """
    # We a assume a single result
    ru_keys = [k for k in pur["protocol_result"]["data"].keys()]
    if len(ru_keys) > 1:
        errmsg = "Too many keys in the Protocol Unit Result file"
        raise ValueError(errmsg)

    units = pur["protocol_result"]["data"][ru_keys[0]]

    if len(units) > 1:
        errmsg = "Too many units found"
        raise ValueError(errmsg)

    if "stateA" not in units[0]["inputs"]:
        return None, None

    stateA = ChemicalSystem.from_dict(units[0]["inputs"]["stateA"])
    stateB = ChemicalSystem.from_dict(units[0]["inputs"]["stateB"])
    ligmap = LigandAtomMapping.from_dict(units[0]["inputs"]["ligandmapping"])

    if any([isinstance(comp, gufe.ProteinComponent) for comp in stateA.components.values()]):
        phase = "complex"
    else:
        phase = "solvent"

    name = f"{phase}_{ligmap.componentA.name}_{ligmap.componentB.name}"

    transform = Transformation(
        stateA=stateA,
        stateB=stateB,
        mapping=ligmap,
        name=name,
        protocol=None,
    )

    return transform, phase


def _check_and_deduplicate_transforms(
    transforms_dict: dict[str, list[Transformation]]
):
    """
    Traverse through a dictionary of transformations keyed
    by the transformation name.

    Check:
      * There are a total of 3 transformations per entry
      * The 3 transformations are the same

    Returns
    -------
    alchemical_network : gufe.AlchemicalNetwork
      A deduplicated alchemical network.
    """
    transform_list = []

    for t_name, t_list in transforms_dict.items():
        if len(t_list) != 3:
            # TODO: we can turn this into a warning message if it's too painful
            errmsg = (
                f"Too few transformations found for {t_name} "
                "this indicates a partially completed set of results. "
                "Please ensure that your input network is finished and "
                "any reproducible partial failures have been removed."
            )
            raise ValueError(errmsg)
        if not all(a == t_list[0] for a in t_list):
            errmsg = (
                f"Transformations for {t_name} do not match "
                "this may indicate some kind of corruption of the alchemical "
                "network. Please contact the OpenFE team."
            )
            raise ValueError(errmsg)

        transform_list.append(t_list[0])

    return AlchemicalNetwork(transform_list)


def parse_results(
    result_files: list[str],
    input_ligand_network: LigandNetwork,
) -> AlchemicalNetwork:
    """
    Create an AlchemicalNetwork from a set of input JSON files.
    """
    # All the transforms, triplicated
    all_transforms_dict = {}

    # Fail if we have no inputs
    if len(result_files) == 0:
        errmsg = "No results .json files were provided"
        raise ValueError(errmsg)

    for resfile in result_files:
        ru = _get_check_results_json(resfile)

        if ru is None:
            continue

        # TODO: fix this, for now this will just fail
        transform, phase = get_transformation(ru)

        if transform is None:
            # We delete the inputs data
            transform, phase = get_transformation_alternate(
                ru, input_ligand_network
            )

        if transform.name in all_transforms_dict:
            all_transforms_dict[transform.name].append(transform)
        else:
            all_transforms_dict[transform.name] = [transform]

    alchemical_network = _check_and_deduplicate_transforms(all_transforms_dict)

    return alchemical_network


def alchemical_network_to_ligand_network(alchemical_network) -> LigandNetwork:
    """
    Create a LigandNetwork from an AlchemicalNetwork
    """
    edges = []
    for e in alchemical_network.edges:
        edges.append(e.mapping)
    # Only add edges if mappings are present twice, meaning both solvent and complex
    # phases are present in the AlchemicalNetwork
    edges = [e for e in edges if edges.count(e) == 2]
    network = LigandNetwork(edges=set(edges))
    return network


def decompose_disconnected_ligand_network(
    ligand_network: LigandNetwork,
) -> list[LigandNetwork]:
    """
    Get a list of connected networks from an input
    disconnected network
    """
    # build alchemical sub-networks
    g = ligand_network.graph
    g = g.to_undirected()
    all_edges = list(ligand_network.edges)
    sub_network_nodes = list(nx.connected_components(g))

    ligand_sub_networks = []
    for i, snn in enumerate(sub_network_nodes):
        # collect edges
        sub_edges = []
        for e in all_edges:
            if any(n in (e.componentA, e.componentB) for n in snn):
                sub_edges.append(e)

        ligand_sub_network = LigandNetwork(edges=sub_edges)
        ligand_sub_networks.append(ligand_sub_network)

    return ligand_sub_networks


def get_new_network_connections(
    ligand_sub_networks: list[LigandNetwork],
    input_ligand_network: LigandNetwork,
) -> LigandNetwork:

    mapper = KartografAtomMapper()
    scorer = lomap_scorers.default_lomap_score

    # Create a maximal network
    print("LOG: generating maximal network -- this may take time")
    max_network = openfe.setup.ligand_network_planning.generate_maximal_network(
        input_ligand_network.nodes,
        mapper,
        scorer,
        True
    )

    input_edges = input_ligand_network.edges

    # Adding reverse input edges, in order to avoid these transfomations totally.
    # CAVEAT: the order in the componentA_to_componentB dict will not always match the
    # mapping atom order created elsewhere, therefore we can't use these to
    # check for edge equality, but must go through the components (see below)
    input_edges = input_edges.union(
        set(
            [
                LigandAtomMapping(
                    componentA=m.componentB,
                    componentB=m.componentA,
                    componentA_to_componentB={
                        v: k for k, v in m.componentA_to_componentB.items()
                    },
                )
                for m in list(input_edges)
            ]
        )
    )

    # Prune out prior edges
    # We have to go through the mapping components instead of the edges
    new_edges = list(max_network.edges).copy()
    for edge in max_network.edges:
        molA = edge.componentA
        molB = edge.componentB
        for e in input_edges:
            if molA == e.componentA and molB == e.componentB:
                new_edges.remove(edge)
                break

    # Prune out self subnetwork edges
    usable_edges = new_edges.copy()
    for sn in ligand_sub_networks:
        for edge in new_edges:
            molA = edge.componentA
            molB = edge.componentB
            if molA in sn.nodes and molB in sn.nodes:
                usable_edges.remove(edge)

    # Sort usable_edges by score (lowest to highest)
    idx_scores = [(index, edge.annotations["score"]) for index, edge in enumerate(usable_edges)]
    sorted_idx_scores = sorted(idx_scores, key=lambda x: x[1])
    sorted_usable_edges = [usable_edges[i[0]] for i in sorted_idx_scores]

    # Calculate redundancies: Adding all possible connecting edges to the
    # respective sub ligand networks
    connections_sub_networks = [[] for _ in ligand_sub_networks]
    for edge in sorted_usable_edges:
        for i, net in enumerate(ligand_sub_networks):
            if edge.componentA in net.nodes or edge.componentB in net.nodes:
                connections_sub_networks[i].append(edge)

    # Get the edges to tape the network
    # We work by pruning, so the edges towards the end (better scored) will
    # be more likely to be added
    tape_edges = []
    for edge in sorted_usable_edges:
        for net in connections_sub_networks:
            if edge in net:
                if len(net) <= 2 and edge not in tape_edges:
                    tape_edges.append(edge)
                else:
                    net.remove(edge)

    # Get a list of the old edges (from the broken network)
    old_edges = [list(net.edges) for net in ligand_sub_networks if len(net.edges) > 0]
    old_edges = [item for list in old_edges for item in list]

    # Check if the broken network is connected with the new tapes
    concatenated_edges = old_edges + tape_edges
    concatenated_network = LigandNetwork(
        nodes=input_ligand_network.nodes,
        edges=concatenated_edges,
    )
    if not concatenated_network.is_connected():
        raise ValueError(
            "The attempt to fix the broken network failed!"
            "Please contact the OpenFE team if you experience this."
        )

    return LigandNetwork(edges=tape_edges)


def get_alchemical_charge_difference(mapping) -> int:
    """
    Checks and returns the difference in formal charge between state A and B.

    Parameters
    ----------
    mapping : dict[str, ComponentMapping]
      Dictionary of mappings between transforming components.

    Returns
    -------
    int
      The formal charge difference between states A and B.
      This is defined as sum(charge state A) - sum(charge state B)
    """
    chg_A = Chem.rdmolops.GetFormalCharge(
        mapping.componentA.to_rdkit()
    )
    chg_B = Chem.rdmolops.GetFormalCharge(
        mapping.componentB.to_rdkit()
    )

    return chg_A - chg_B


def get_settings():
    """
    Utility method for getting RFEProtocol settings for non charge changing
    transformations.
    These settings mostly follow defaults but use the newest OpenFF 2.2.
    """
    # Are there additional settings we should specify here?
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Should we use this new OpenFF version or the default?
    settings.forcefield_settings.small_molecule_forcefield = 'openff-2.2.0'
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def get_settings_charge_changes():
    """
    Utility method for getting RFEProtocol settings for charge changing
    transformations.

    These settings mostly follow defaults but use longer
    simulation times, more lambda windows and an alchemical ion.
    """
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Should we use this new OpenFF version or the default?
    settings.forcefield_settings.small_molecule_forcefield = 'openff-2.2.0'
    settings.alchemical_settings.explicit_charge_correction = True
    settings.simulation_settings.production_length = 20 * unit.nanosecond
    settings.simulation_settings.n_replicas = 22
    settings.lambda_settings.lambda_windows = 22
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def get_fixed_alchemical_network(ducktape_network, alchemical_network):
    """
    Create an alchemical networks with only the missing edges.
    """
    solv = openfe.SolventComponent()

    for node in alchemical_network.nodes:
        prot_comps = [isinstance(comp, ProteinComponent) for comp in node.components.values()]
        if len(prot_comps) > 0:
            prot = prot_comps[0]
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
        charge_difference = get_alchemical_charge_difference(mapping)
        if abs(charge_difference) > 1e-3:
            # Raise a warning that a charge changing transformation is included
            # in the network
            wmsg = (
                "Charge changing transformation between ligands "
                f"{mapping.componentA.name} and "
                f"{mapping.componentB.name}. "
                "A more expensive protocol with 22 lambda windows, "
                "sampled for 20 ns each, will be used here."
            )
            warnings.warn(wmsg)
            # Get settings for charge changing transformations
            rfe_settings = get_settings_charge_changes()
        else:
            rfe_settings = get_settings()
        for leg in ["solvent", "complex"]:
            # use the solvent and protein created above
            sysA_dict = {"ligand": mapping.componentA, "solvent": solv}
            sysB_dict = {"ligand": mapping.componentB, "solvent": solv}

            if leg == "complex":
                sysA_dict["protein"] = prot
                sysB_dict["protein"] = prot
                if len(cofactors) > 0:

                    for cofactor, entry in zip(cofactors_smc, string.ascii_lowercase):
                        cofactor_name = f"cofactor_{entry}"
                        sysA_dict[cofactor_name] = cofactor
                        sysB_dict[cofactor_name] = cofactor

            sysA = openfe.ChemicalSystem(sysA_dict)
            sysB = openfe.ChemicalSystem(sysB_dict)

            name = f"{leg}_{mapping.componentA.name}_" f"{mapping.componentB.name}"

            rbfe_protocol = RelativeHybridTopologyProtocol(settings=rfe_settings)
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


def fix_network(
    result_files: list[str],
    input_alchem_network_file: pathlib.Path,
    output_alchemical_network_folder: pathlib.Path,
):
    print("Parsing input files:")
    # Parse Ligand Network
    print("LOG: Reading input Alchemical Network")
    input_alchem_network = parse_alchemical_network(input_alchem_network_file)
    input_ligand_network = alchemical_network_to_ligand_network(input_alchem_network)

    # Parse Alchemical Network
    print("LOG: Reading alchemical network result files")
    res_alchemical_network = parse_results(result_files, input_ligand_network)
    res_ligand_network = alchemical_network_to_ligand_network(
        res_alchemical_network
    )

    # decompose disconnected sets
    ligand_sub_networks = decompose_disconnected_ligand_network(
        res_ligand_network
    )

    # get the ligands which are missing
    missing_ligands = list(
        input_ligand_network.nodes.difference(res_ligand_network.nodes)
    )

    # turn the missing ligands into spare ligand networks
    for missing_ligand in missing_ligands:
        ligand_sub_networks.append(LigandNetwork(nodes=[missing_ligand], edges=[]))

    # Sort the sub-networks by their size in order of biggest to smallest
    ligand_sub_networks = list(
        sorted(ligand_sub_networks, key=lambda sn: len(sn.nodes), reverse=True)
    )

    # check status
    print("LOG: Reporting input file contents: ")
    print(f"\tPlanned input  no. ligands: {len(input_ligand_network.nodes)}\n")
    print(f"\tPlanned input  no. connections: {len(input_ligand_network.edges)}\n")
    print(f"\tSimulation results no. ligands: {len(res_ligand_network.nodes)}\n")
    print(f"\tSimulation results no. connections: {len(res_ligand_network.edges)}\n")
    print(f"\tMissing ligands in simulation results: {len(missing_ligands)}\n")
    #ToDo: If there is no disconnected network this would return 1 which is confusing
    print(f"\tDisconnected networks which need patching: {len(ligand_sub_networks)}\n")
    print(
        f"\tLigands in each disconnected network: {[len(sn.nodes) for sn in ligand_sub_networks]}"
    )
    print("\n")

    if len(ligand_sub_networks) == 1:
        # Return gracefully, don't fail!
        message = ("Did not find disconnected components in alchemical "
                   "network, nothing to do here!")
        print(message)
        return None

    print()

    # Tape the networks together
    print("LOG: finding additional connections to merge the broken networks:")
    network_connections = get_new_network_connections(
        ligand_sub_networks=ligand_sub_networks,
        input_ligand_network=input_ligand_network,
    )

    print(f"\tLOG: Generated new ligand transform edges: {len(network_tapes.edges)}\n")

    # write out
    print("LOG: writing out network of additional edges to fix the network:")
    output_alchemical_network_folder.mkdir(exist_ok=False, parents=True)
    taped_alchemical_network = get_fixed_alchemical_network(
        network_connections, input_alchem_network
    )
    alchemical_network_json_fp = (
        output_alchemical_network_folder / "alchemical_network.json"
    )
    json.dump(
        taped_alchemical_network.to_dict(),
        alchemical_network_json_fp.open(mode="w"),
        cls=JSON_HANDLER.encoder,
    )

    # Write out each transformation
    # Create a subdirectory for the transformations
    transforms_dir = pathlib.Path(output_alchemical_network_folder / "transformations")
    transforms_dir.mkdir(exist_ok=True, parents=True)

    for transform in taped_alchemical_network.edges:
        transform.dump(transforms_dir / f"{transform.name}.json")

    print("Done!")


def cli_fix_network():
    import argparse
    parser = argparse.ArgumentParser(
        description="Fix broken alchemical network"
    )
    parser.add_argument(
        "--input_alchem_network_file", type=pathlib.Path, help="Path to the input alchemical network",
    )
    parser.add_argument(
        "--output_extra_transformations", type=pathlib.Path, help="Path to where we will write out extra transformations to fix the network",
    )
    parser.add_argument(
        "--result_files", nargs="+", help="Results JSON file(s)"
    )


    args = parser.parse_args()
    fix_network(
        result_files=args.result_files,
        input_alchem_network_file=args.input_alchem_network_file,
        output_alchemical_network_folder=args.output_extra_transformations,
    )


# main Runs
if __name__ == "__main__":
    cli_fix_network()