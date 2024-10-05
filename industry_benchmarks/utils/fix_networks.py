import json
import click
import pathlib
from typing import Iterable
import numpy as np
import itertools
from glob import glob
import networkx as nx
import warnings

import gufe
import pytest
from gufe.tokenization import JSON_HANDLER
import openfe
from openfe import ChemicalSystem, LigandAtomMapping, Protocol
from openfe import Transformation, AlchemicalNetwork, LigandNetwork
from openfe import SolventComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe.setup import KartografAtomMapper, lomap_scorers
from . import plan_rbfe_network


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
    if not os.path.exists(json_file):
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


def parse_results(
    result_files: list[str],
    input_ligand_network: LigandNetwork,
) -> AlchemicalNetwork:
    """
    Create an AlchemicalNetwork from a set of input JSON files.
    """
    nodes = []
    edges = {}
    edge_count = {}

    # Fail if we have no inputs
    if len(result_files) == 0:
        errmsg = "No results .json files were provided"
        raise ValueError(errmsg)

    for resfile in result_files:
        ru = _get_check_results_json(resfile)

        if ru is None:
            continue

        for k in res["protocol_result"]["data"]:
            units = res["protocol_result"]["data"][k]

            name = units[0]["name"].split(" repeat")[0]
            # TODO: use protocol_result instead - that's not cleaned up
            stateA = ChemicalSystem.from_dict(units[0]["inputs"]["stateA"])
            stateB = ChemicalSystem.from_dict(units[0]["inputs"]["stateB"])
            ligmap = LigandAtomMapping.from_dict(
                units[0]["inputs"]["ligandmapping"]
            )

            # Count if both transformations present.
            if name not in edge_count:
                edge_count[name] = {"protein": 0, "solvent": 0}

            if "protein" in stateA.components or "solvent" not in stateA.components:
                edge_count[name]["protein"] += 1
            else:
                edge_count[name]["solvent"] += 1

            t = Transformation(
                stateA=stateA,
                stateB=stateB,
                mapping=ligmap,
                name=name,
                protocol=None,
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
    """
    Create a LigandNetwork from an AlchemicalNetwork
    """
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


#def ligand_network_to_networkx(ligand_network) -> nx.Graph:
#
#    edges = ligand_network.edges
#    wedges = []
#    for edge in list(edges):
#        wedges.append([edge.componentA, edge.componentB, edge.annotations["score"]])
#
#    g = nx.Graph(nodes=ligand_network.nodes)
#    g.add_weighted_edges_from(ebunch_to_add=wedges)
#    return g


def decompose_disconnected_ligand_network(
    ligand_network: LigandNetwork,
) -> list[LigandNetwork]:
    """
    Get a list of connected networks from an input
    disconnected network
    """
    # build alchemical sub-networks
    g = ligand_network.graph
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
) -> LigandNetwork:

    mapper = KartografAtomMapper()
    scorer = lomap_scorers.default_lomap_score

    # Create a maximal networka
    print("LOG: generating maximal network -- this may take time")
    max_network = openfe.setup.generate_maximal_network(
        input_ligand_network.nodes,
        mapper,
        scorer,
        True
    )

    in_edges = input_ligand_network.edges

    # Adding reverse input edges, in order to avoid these transfomations totally.
    # CAVE: the order in the componentA_to_componentB dict will not always match the
    # mapping atom order created elsewhere, therefore we can't use these to
    # check for edge equality, but must go through the components (see below)
    in_edges = in_edges.union(
        set(
            [
                LigandAtomMapping(
                    componentA=m.componentB,
                    componentB=m.componentA,
                    componentA_to_componentB={
                        v: k for k, v in m.componentA_to_componentB.items()
                    },
                )
                for m in list(in_edges)
            ]
        )
    )

    # Prune out prior edges
    # We have to go through the mapping components instead of the edges
    new_edges = list(max_network.edges).copy()
    for edge in max_network.edges:
        molA = edge.componentA
        molB = edge.componentB
        for e in in_edges:
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
    scores = [edge.annotations["score"] for edge in usable_edges]
    sorted_score_indices = np.argsort(scores)
    sorted_usable_edges = [usable_edges[i] for i in sorted_score_indices]

    # Calculate redundancies: Adding all possible connecting edges to the
    # respective sub ligand networks
    connections_sub_networks = [[] for _ in ligand_sub_networks]
    for edge in sorted_usable_edges:
        for i, net in enumerate(ligand_sub_networks):
            if edge.componentA in net.nodes or edge.componentB in net.nodes:
                connections_sub_networks[i].append(edge)

    # Get the edges to tape the network
    tape_edges = []
    for edge in sorted_usable_edges:
        for i, net in enumerate(connections_sub_networks):
            if edge in net:
                if len(net) <= 2 and edge not in tape_edges:
                    tape_edges.append(edge)
                else:
                    connections_sub_networks[i].remove(edge)

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

    tape_nodes = set([n for e in tape_edges for n in [e.componentA, e.componentB]])

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
            wmsg = (
                "Charge changing transformation between ligands "
                f"{mapping.componentA.name} and "
                f"{mapping.componentB.name}. "
                "A more expensive protocol with 22 lambda windows, "
                "sampled "
                "for 20 ns each, will be used here."
            )
            warnings.warn(wmsg)
            # Get settings for charge changing transformations
            rfe_settings = plan_rbfe_network.get_settings_charge_changes()
        else:
            rfe_settings = plan_rbfe_network.get_settings()
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


def do_taping(
    result_files: list[str],
    input_alchem_network_file: pathlib.Path,
    output_alchemical_network_folder: pathlib.Path,
):
    print("Parsing input files:")
    # Parse Ligand Network
    print("Log: Reading input Alchemical Network")
    input_alchem_network = parse_alchemical_network(input_alchem_network_file)
    input_ligand_network = alchemical_network_to_ligand_network(input_alchem_network)

    # Parse Alchemical Network
    print("Log: Reading result files")
    res_alchemical_network = parse_results(
        result_files_regex, input_ligand_network
    )
    res_ligand_network = alchemical_network_to_ligand_network(res_alchemical_network)

    if res_ligand_network.is_connected:
        # TODO: 
        print("LOG: results network is connected, nothing to do here")
        return None

    # decompose disconnected sets
    ligand_sub_networks = decompose_disconnected_ligand_network(res_ligand_network)

    # Add missing nodes to ligand_sub_networks if ligand network provided
    # get network node differences
    alchemical_network_ligands = set([n for n in res_ligand_network.nodes])
    missing_commponents = list(
        input_ligand_network.nodes.difference(res_ligand_network.nodes)
    )

    for missing_commponent in missing_commponents:
        ligand_sub_networks.append(LigandNetwork(nodes=[missing_commponent], edges=[]))
    ligand_sub_networks = list(
        sorted(ligand_sub_networks, key=lambda sn: len(sn.nodes), reverse=True)
    )

    # check status
    print(f"Inspecting input files: ")
    print(f"\tinput LigandNetwork  nodes: {len(input_ligand_network.nodes)}")
    print(f"\tinput LigandNetwork  edges: {len(input_ligand_network.edges)}")
    print(f"\tresult LigandNetwork nodes: {len(res_ligand_network.nodes)}")
    print(f"\tresult LigandNetwork edges: {len(res_ligand_network.edges)}")
    print()

    print(f"Found the following problems: ")
    print(f"\tmissing components: {len(missing_commponents)}")
    print(f"\tdisconnected networks: {len(ligand_sub_networks)}")
    print(
        f"\tdisconnected networks sizes: {[len(sn.nodes) for sn in ligand_sub_networks]}"
    )

    if len(ligand_sub_networks) == 1:
        raise ValueError("Did not find disconnected components in alchemical network!")

    print()

    # Tape the networks together
    print("Taping the Networks together:")
    network_tapes = get_new_network_tapes(
        ligand_sub_networks=ligand_sub_networks,
        input_ligand_network=input_ligand_network,
    )

    print(f"\tGenerated new taping LigandMapping edges: {len(network_tapes.edges)}")
    print(f"\tUsing ligand nodes for taping: {len(network_tapes.nodes)}")
    print()

    # write out
    print("Write out the tapes for the network:")
    output_alchemical_network_folder.mkdir(exist_ok=False, parents=True)
    taped_alchemical_network = get_taped_alchemical_network(
        network_tapes, input_alchem_network
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


def cli_do_taping():
    import argparse
    parser = argparse.ArgumentParser(
        description="Fix broken alchemical network"
    )
    parser.add_argument(
        "result_files", metavar="JSON_FILE", nargs="+", help="Results JSON file(s)"
    )
    parser.add_argument(
        "input_alchem_network_file", type=pathlib.Path, help="Path to the input alchemical network",
    )
    parser.add_argument(
        "output_extra_transformations", type=Pathlib.Path, help="Path to where we will write out extra transformations to fix the network",
    )

    args = parser.parse_args()
    do_taping(
        result_files=args.result_files,
        input_alchem_network_file=args.input_alchem_network_file,
        output_alchemical_network_folder=args.output_extra_transformations,
    )


# main Runs
if __name__ == "__main__":
    cli_do_taping()
