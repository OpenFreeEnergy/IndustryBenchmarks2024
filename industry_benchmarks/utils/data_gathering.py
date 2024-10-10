import click
import pathlib
import json
from rdkit import Chem
import gufe
import openfe
from openfe import LigandNetwork
from kartograf.atom_mapping_scorer import (
    MappingRMSDScorer, MappingShapeOverlapScorer, MappingVolumeRatioScorer,
)


def parse_ligand_network(
    input_ligand_network_path: str,
) -> LigandNetwork:
    """
    Read a LigandNetwork
    """
    with open(input_ligand_network_path) as f:
        graphml = f.read()

    network = gufe.LigandNetwork.from_graphml(graphml)
    return network


def get_blinded_transformation_network(
    input_ligand_network: LigandNetwork,
):
    """
    Blinded transformation network: Names of nodes and how they are connected.
    Remove smcs, mappings
    """

    return


def get_lomap_score(mapping):
    score = openfe.setup.lomap_scorers.default_lomap_score(mapping)
    return score


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


def get_shape_score(mapping):
    shape_scorer = MappingShapeOverlapScorer()
    score = shape_scorer(mapping=mapping)
    return score


def get_volume_score(mapping):
    volume_scorer = MappingVolumeRatioScorer()
    score = volume_scorer(mapping=mapping)
    return score


def get_mapping_RMSD_score(mapping):
    rmsd_scorer = MappingRMSDScorer()
    score = rmsd_scorer(mapping=mapping)
    return score


def get_number_heavy_dummy_heavy_core_atoms():
    return


def get_fingerprint_similarity_score():
    """
    2/3D fingerprint similarity scores (e.g. Tanimoto)
    """
    return


def get_changing_number_rotatable_bonds():
    """
    Number of rotatable bonds in the non-core region
    """
    return


def get_changing_number_rings():
    """
    Number of ring systems in the non-core region
    """
    return


def gather_transformation_scores(
    input_ligand_network: LigandNetwork,
) -> dict[str, dict[str, int]]:
    """
    Gather all scores for the transformations in a dict.
    Transformations information
        Lomap score
        formal charge score/number of charge changes
        Shape score
        Volume score
        RMSD (mapping RMSD score)
        Number of heavy dummy and core atoms
        2/3D fingerprint similarity scores (e.g. Tanimoto)
        changing number of rotatable bonds
        changing number of rings
    Returns
    -------
    dict[str, dict[str, int]]
      Dictionary of the edge name and a dictionary of score name and value.
    """
    transformations_scores = {}
    for edge in input_ligand_network.edges:

        name = f'edge_{edge.componentA.name}_{edge.componentB.name}'
        transformations_scores[name] = {}
        edge_scores = {}
        lomap_score = get_lomap_score(edge)
        edge_scores["lomap_score"] = lomap_score
        alchemical_charge_difference = get_alchemical_charge_difference(edge)
        edge_scores["alchemical_charge_difference"] = alchemical_charge_difference
        shape_score = get_shape_score(edge)
        edge_scores["shape_score"] = shape_score
        volume_score = get_volume_score(edge)
        edge_scores["volume_score"] = volume_score
        mapping_rmsd_score = get_mapping_RMSD_score(edge)
        edge_scores["mapping_rmsd_score"] = mapping_rmsd_score

        transformations_scores[name] = edge_scores

    return transformations_scores


def get_number_rotatable_bonds(smc):
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(m, strict=True)
    return num_rotatable_bonds

def get_number_ring_systems(smc):
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_rings = Chem.rdMolDescriptors.CalcNumRings(m)
    return num_rings

def get_number_heavy_atoms(smc):
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_heavy_atoms = Chem.rdMolDescriptors.CalcNumHeavyAtoms(m)
    return num_heavy_atoms

def get_system_element_count(smc):
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    atomic_numbers = [atom.GetAtomicNum() for atom in m.GetAtoms()]
    return len(set(atomic_numbers))

def gather_ligand_scores(
    input_ligand_network: LigandNetwork,
) -> dict[str, dict[str, int]]:
    """
    Gather all scores for the ligands in a dict.
    Ligand information
        Number of rotatable bonds
        Number of rings
        Number of heavy atoms
        system element counts
    Returns
    -------
    dict[str, dict[str, int]]
      Dictionary of the ligand name and a dictionary of score name and value.
    """
    all_ligand_scores = {}
    for node in input_ligand_network.nodes:
        name = f'ligand_{node.name}'
        all_ligand_scores[name] = {}
        ligand_scores = {}
        num_rotatable_bonds = get_number_rotatable_bonds(node)
        ligand_scores["num_rotatable_bonds"] = num_rotatable_bonds
        num_rings = get_number_ring_systems(node)
        ligand_scores["num_rings"] = num_rings
        num_heavy_atoms = get_number_heavy_atoms(node)
        ligand_scores["num_heavy_atoms"] = num_heavy_atoms
        num_elements = get_system_element_count(node)
        ligand_scores["num_elements"] = num_elements

        all_ligand_scores[name] = ligand_scores

    return all_ligand_scores


@click.command
@click.option(
    '--input_ligand_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./alchemicalNetwork/ligand_network.graphml"),
    required=True,
    help=("Path to the ligand_network.graphml file that was used to run these "
         "simulations."),
)
@click.option(
    '--output_scores',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./scores.json"),
    required=True,
    help=("Path to the JSON file that stores all the scores and metrics for "
          "this dataset."),
)
@click.option(
    '--fixed_ligand_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./ligand_network.graphml"),
    required=False,
    help=("Only needed when a broken network was fixed with additional edges. "
          "Path to the ligand_network.graphml file that was used to run the "
         "simulations of fixing the network."),
)
def gather_data(
    input_ligand_network,
    output_scores,
    fixed_ligand_network,
):
    """
    Function that gathers all the data.
    """
    ligand_network = parse_ligand_network(input_ligand_network)
    transformation_scores = gather_transformation_scores(ligand_network)
    ligand_scores = gather_ligand_scores(ligand_network)
    # Create a single dict of all scores
    scores = {
        "transformation_scores": transformation_scores,
        "ligand_scores": ligand_scores,
    }
    # Save this to json
    file = pathlib.Path(output_scores)
    with open(file, mode='w') as f:
        json.dump(scores, f)


if __name__ == "__main__":
    gather_data()
