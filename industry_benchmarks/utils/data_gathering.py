import gufe
from openfe import LigandNetwork


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


def get_lomap_score():
    return


def get_formal_charge_score():
    """
    Whether the transformation undergoes a change in net charge.
    """
    return


def get_shape_score():
    return


def get_volume_score():
    return


def get_mapping_RMSD_score():
    return

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

    return


def get_number_rotatable_bonds():
    return

def get_number_ring_systems():
    return

def get_number_heavy_atoms():
    return

def get_system_element_count():
    return

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
    return


def gather_data(
    input_ligand_network_path: str,
):
    """
    Function that gathers all the data.
    """
    ligand_network = parse_ligand_network(input_ligand_network_path)
    transformation_scores = gather_transformation_scores(ligand_network)
    # Save this to json

    ligand_scores = gather_ligand_scores(ligand_network)
    # Save this to json
