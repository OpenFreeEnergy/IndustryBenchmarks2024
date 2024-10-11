from ..data_gathering import (
    get_lomap_score,
    get_blinded_transformation_network,
    parse_ligand_network,
    get_shape_score,
    get_volume_score,
    get_mapping_RMSD_score,
    get_number_heavy_dummy_heavy_core_atoms,
    gather_transformation_scores,
    get_number_rotatable_bonds,
    get_number_ring_systems,
    get_number_heavy_atoms,
    get_system_element_count,
    gather_ligand_scores,
    get_changing_number_rings,
    get_changing_number_rotatable_bonds
)
import pytest
from importlib import resources
from gufe import LigandNetwork
import numpy as np

@pytest.fixture
def cmet_ligand_network() -> LigandNetwork:
    with resources.files("utils.tests.data") as d:
        yield parse_ligand_network(str(d / "cmet_results/alchemicalNetwork/ligand_network.graphml"))

# test all metrics independently
def test_blinded_network(cmet_ligand_network):
    """Test extracting blinded data from the ligand network, this is just the names of nodes and the edges."""
    transform_info = get_blinded_transformation_network(cmet_ligand_network)
    for node in cmet_ligand_network.nodes:
        assert node.name in transform_info["nodes"]
    for edge in cmet_ligand_network.edges:
        assert (edge.componentA.name, edge.componentB.name) in transform_info["edge"]


def test_get_shape_score(cmet_ligand_network):
    expected_scores = [0.16995560327261203, 0.18790420591686238, 0.21368186150398924, 0.24393553760289846]
    scores = sorted([get_shape_score(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(scores, expected_scores)


def test_get_volume_score(cmet_ligand_network):
    expected_scores = [0.4630927 , 0.46839382, 0.49807746, 0.75730246]
    scores = sorted([get_volume_score(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(scores, expected_scores)

def test_mapping_rmsd_score(cmet_ligand_network):
    expected_scores = [0.55890129, 0.61935956, 0.6212269 , 0.62804141]
    scores = sorted([get_mapping_RMSD_score(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(scores, expected_scores)

def test_number_heavy_dummy_atoms(cmet_ligand_network):
    expected_nos = [[19, 9, 8], [22, 6, 6], [22, 6, 7], [28, 1, 9]]
    atom_counts = sorted([list(get_number_heavy_dummy_heavy_core_atoms(edge)) for edge in cmet_ligand_network.edges])
    assert np.allclose(atom_counts, expected_nos)

def test_get_lomap_score(cmet_ligand_network):
    expected_scores = [0.10025884, 0.16529889, 0.18268352, 0.2865048 ]
    scores = sorted([get_lomap_score(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(scores, expected_scores)

def test_gather_transfer_info(cmet_ligand_network):
    """Test gathering all of the scores for a network"""
    transformation_scores = gather_transformation_scores(cmet_ligand_network)
    # make sure all the edges are present
    assert len(transformation_scores) == 4
    # for each edge make sure that all the expected data is present
    for edge_info in transformation_scores.values():
        assert "lomap_score" in edge_info
        assert "alchemical_charge_difference" in edge_info
        assert "shape_score" in edge_info
        assert "volume_score" in edge_info
        assert "mapping_rmsd_score" in edge_info
        assert "num_heavy_core" in edge_info
        assert "num_heavy_dummy_A" in edge_info
        assert "num_heavy_dummy_B" in edge_info
        assert "difference_num_rings_AB" in edge_info
        assert "difference_num_rot_bonds_AB" in edge_info


def test_number_of_rotor_bonds(cmet_ligand_network):
    expected_rotors = [5, 5, 5, 6, 8]
    rotors = sorted([get_number_rotatable_bonds(node) for node in cmet_ligand_network.nodes])
    assert np.allclose(rotors, expected_rotors)

def test_number_of_rings(cmet_ligand_network):
    expected_rings = [3, 3, 4, 4, 5]
    rings = sorted([get_number_ring_systems(node) for node in cmet_ligand_network.nodes])
    assert np.allclose(rings, expected_rings)

def test_number_heavy_atoms(cmet_ligand_network):
    expected_heavy_atoms = [27, 28, 28, 29, 37]
    heavy_atoms = sorted([get_number_heavy_atoms(node) for node in cmet_ligand_network.nodes])
    assert np.allclose(heavy_atoms, expected_heavy_atoms)

def test_number_of_elements(cmet_ligand_network):
    expected_element_numbers = [4, 5, 5, 5, 5]
    element_numbers = sorted([get_system_element_count(node) for node in cmet_ligand_network.nodes])
    assert np.allclose(element_numbers, expected_element_numbers)

def test_gather_ligand_scores(cmet_ligand_network):
    all_ligand_scores = gather_ligand_scores(cmet_ligand_network)
    assert len(all_ligand_scores) == 5
    for ligand_data in all_ligand_scores.values():
        assert "num_rotatable_bonds" in ligand_data
        assert "num_rings" in ligand_data
        assert "num_heavy_atoms" in ligand_data
        assert "num_elements" in ligand_data


def test_changing_number_of_rings(cmet_ligand_network):
    expected_ring_changes = [0, 1, 1, 1]
    ring_changes = sorted([get_changing_number_rings(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(ring_changes, expected_ring_changes)

def test_changing_number_of_rotatable_bonds(cmet_ligand_network):
    expected_rotor_changes = [1, 1, 1, 3]
    rotor_changes = sorted([get_changing_number_rotatable_bonds(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(rotor_changes, expected_rotor_changes)