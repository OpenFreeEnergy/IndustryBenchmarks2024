import click
import pathlib
import json
from rdkit import Chem
from rdkit.Chem import AllChem
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
    blinded_network = {}
    blinded_network["nodes"] = [node.name for node in input_ligand_network.nodes]
    blinded_network["edge"] = [(edge.componentA.name, edge.componentB.name) for edge in input_ligand_network.edges]

    return blinded_network


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


def get_number_heavy_dummy_heavy_core_atoms(mapping):
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    mapped_atomIDs = mapping.componentA_to_componentB

    heavy_core_A = []
    heavy_dummy_A = []
    heavy_core_B = []
    heavy_dummy_B = []

    for atom in molA.GetAtoms():
        # For heavy core: Heavy atom indices that are in the mapping:
        if atom.GetIdx() in mapped_atomIDs.keys() and atom.GetAtomicNum() != 1:
            # Check if the atom is mapped to a hydrogen atom
            index_mapping = list(mapped_atomIDs.keys()).index(atom.GetIdx())
            index_atomB = list(mapped_atomIDs.values())[index_mapping]
            atomic_number_mapped_atom_B = molB.GetAtoms()[index_atomB].GetAtomicNum()
            # If the mapped atom in ligand B is a hydrogen, these atoms are unmapped.
            # Therefore we will could the heavy atom as a dummy atom, not core
            if atomic_number_mapped_atom_B == 1:
                heavy_dummy_A.append(atom.GetIdx())
            else:
                heavy_core_A.append(atom.GetIdx())
        # For heavy dummy: Heavy atom indices that are not in the mapping:
        if atom.GetIdx() not in mapped_atomIDs.keys() and atom.GetAtomicNum() != 1:
            heavy_dummy_A.append(atom.GetIdx())

    for atom in molB.GetAtoms():
        # For heavy core: Heavy atom indices that are in the mapping:
        if atom.GetIdx() in mapped_atomIDs.values() and atom.GetAtomicNum() != 1:
            # Check if the atom is mapped to a hydrogen atom
            index_mapping = list(mapped_atomIDs.values()).index(atom.GetIdx())
            index_atomA = list(mapped_atomIDs.keys())[index_mapping]
            atomic_number_mapped_atom_A = molA.GetAtoms()[index_atomA].GetAtomicNum()
            # If the mapped atom in ligand B is a hydrogen, these atoms are unmapped.
            # Therefore we will could the heavy atom as a dummy atom, not core
            if atomic_number_mapped_atom_A == 1:
                heavy_dummy_B.append(atom.GetIdx())
            else:
                heavy_core_B.append(atom.GetIdx())
        # For heavy dummy: Heavy atom indices that are not in the mapping:
        if atom.GetIdx() not in mapped_atomIDs.values() and atom.GetAtomicNum() != 1:
            heavy_dummy_B.append(atom.GetIdx())

    assert len(heavy_core_A) == len(heavy_core_B)

    return len(heavy_core_A), len(heavy_dummy_A), len(heavy_dummy_B)


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
        num_heavy_core, num_heavy_dummy_A, num_heavy_dummy_B = get_number_heavy_dummy_heavy_core_atoms(edge)
        edge_scores["num_heavy_core"] = num_heavy_core
        edge_scores["num_heavy_dummy_A"] = num_heavy_dummy_A
        edge_scores["num_heavy_dummy_B"] = num_heavy_dummy_B
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
    '--output_dir',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("./"),
    required=True,
    help="Path to the output directory that stores all data.",
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
    output_dir,
    fixed_ligand_network,
):
    """
    Function that gathers all the data.
    """
    # Make folder for outputs
    output_dir.mkdir(exist_ok=False, parents=True)

    ligand_network = parse_ligand_network(input_ligand_network)
    transformation_scores = gather_transformation_scores(ligand_network)
    ligand_scores = gather_ligand_scores(ligand_network)
    # Create a single dict of all scores
    scores = {
        "transformation_scores": transformation_scores,
        "ligand_scores": ligand_scores,
    }
    # Save this to json
    file = pathlib.Path(output_dir / 'scores.json')
    with open(file, mode='w') as f:
        json.dump(scores, f)

    blinded_network = get_blinded_transformation_network(ligand_network)
    # Save this to json
    file = pathlib.Path(output_dir / 'blinded_network.json')
    with open(file, mode='w') as f:
        json.dump(blinded_network, f)


if __name__ == "__main__":
    gather_data()
