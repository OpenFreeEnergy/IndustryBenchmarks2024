import json
import pathlib

from ..data_gathering import (
    get_lomap_score,
    get_transformation_network_map,
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
    get_changing_number_rotatable_bonds,
    get_fingerprint_similarity_score,
    get_difference_solvent_accessible_surface_area,
    get_charge_score,
    get_alchemical_charge_difference,
    get_formal_charge,
    gather_data,
    parse_alchemical_network,
    extract_ligand_network,
    load_results_file,
    get_estimate,
    get_transform_name,
    check_network_is_connected,
)
import pytest
from importlib import resources
from gufe import LigandNetwork, AlchemicalNetwork
from gufe.tokenization import JSON_HANDLER
import numpy as np
from click.testing import CliRunner


@pytest.fixture
def cmet_ligand_network() -> LigandNetwork:
    with resources.files("utils.tests.data") as d:
        yield parse_ligand_network(str(d / "cmet_results/alchemicalNetwork/ligand_network.graphml"))

@pytest.fixture
def cmet_results():
    with resources.files("utils.tests.data.cmet_results") as d:
        yield d

@pytest.fixture
def cmet_network():
    with resources.files("utils.tests.data.cmet_results") as d:
        yield str(d / "alchemicalNetwork/alchemical_network.json")

@pytest.fixture
def bace_network():
    with resources.files("utils.tests.data.bace_results") as d:
        yield str(d / "alchemicalNetwork/alchemical_network.json")

@pytest.fixture
def bace_complete_results():
    with resources.files("utils.tests.data.bace_results") as d:
        yield d

@pytest.fixture
def bace_cleaned_result():
    with resources.files("utils.tests.data.bace_results") as d:
        yield str(d / "cleaned_results/complex_spiro10_spiro6.json")

@pytest.fixture
def bace_full_results():
    with resources.files("utils.tests.data.bace_full_results") as d:
        yield d


# test all metrics independently
def test_blinded_network(cmet_ligand_network):
    """Test extracting blinded data from the ligand network, this is just the names of nodes and the edges."""
    transform_info = get_transformation_network_map(cmet_ligand_network)
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

def test_get_charge_score(cmet_ligand_network):
    expected_scores = [1, 1, 1, 1]
    scores = sorted([get_charge_score(edge) for edge in cmet_ligand_network.edges])
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
        assert "morgan_tanimoto_similarity" in edge_info
        assert "difference_solvent_accessible_surface_area" in edge_info
        assert "charge_score" in edge_info


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

def test_formal_charge(cmet_ligand_network):
    expected_charges = [0, 0, 0, 0, 0]
    charges = sorted([get_formal_charge(node) for node in cmet_ligand_network.nodes])
    assert np.allclose(charges, expected_charges)

def test_gather_ligand_scores(cmet_ligand_network):
    all_ligand_scores = gather_ligand_scores(cmet_ligand_network)
    assert len(all_ligand_scores) == 5
    for ligand_data in all_ligand_scores.values():
        assert "num_rotatable_bonds" in ligand_data
        assert "num_rings" in ligand_data
        assert "num_heavy_atoms" in ligand_data
        assert "num_elements" in ligand_data
        assert "solvent_accessible_surface_area" in ligand_data
        assert "formal_charge" in ligand_data


def test_changing_number_of_rings(cmet_ligand_network):
    expected_ring_changes = [0, 1, 1, 1]
    ring_changes = sorted([get_changing_number_rings(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(ring_changes, expected_ring_changes)

def test_changing_number_of_rotatable_bonds(cmet_ligand_network):
    expected_rotor_changes = [1, 1, 1, 3]
    rotor_changes = sorted([get_changing_number_rotatable_bonds(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(rotor_changes, expected_rotor_changes)


def test_fingerprint_similarity(cmet_ligand_network):
    expected_similarity = [0.45217391, 0.48543689, 0.50505051, 0.81927711]
    similarity = sorted([get_fingerprint_similarity_score(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(similarity, expected_similarity)


def test_get_sasa_diff(cmet_ligand_network):
    expected_sasa_diff = [1.3761717, 14.20480967, 28.02287447, 148.0369631]
    sasa_diff = sorted([get_difference_solvent_accessible_surface_area(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(sasa_diff, expected_sasa_diff)


def test_charge_diff_no_diff(cmet_ligand_network):
    expected_diff = [0, 0, 0, 0]
    charge_diff = sorted([get_alchemical_charge_difference(edge) for edge in cmet_ligand_network.edges])
    assert np.allclose(charge_diff, expected_diff)


def test_parse_alchemical_network(bace_network):
    network = parse_alchemical_network(bace_network)
    assert len(network.nodes) == 18
    assert len(network.edges) == 22
    assert isinstance(network, AlchemicalNetwork)


def test_extract_ligand_network(bace_network):
    network = parse_alchemical_network(bace_network)
    ligand_network = extract_ligand_network(network)
    assert isinstance(ligand_network, LigandNetwork)
    # make sure they have the correct nodes and edges
    assert len(network.nodes) == len(ligand_network.nodes) * 2
    assert len(network.edges) == len(ligand_network.edges) * 2


def test_load_results_file_wrong_file(bace_network, capsys):
    """Make sure None is returned if we try and load the wrong file."""
    result = load_results_file(bace_network)
    assert result is None
    assert "is a AlchemicalNetwork json, skipping" in capsys.readouterr().out

def test_load_results_not_clean(bace_complete_results):
    """Make sure an error is raised if the results have not been cleaned up."""
    with pytest.raises(ValueError, match="has not been cleaned up please make sure you run the `results_cleanup.py` script first"):
        _ = load_results_file(str(bace_complete_results / "results_0" / "solvent_spiro2_spiro1.json"))


def test_get_estimate(cmet_results):
    result = load_results_file(str(cmet_results / "results_0" / "lig_CHEMBL3402745_200_5_solvent_lig_CHEMBL3402744_300_4_solvent_solvent.json"))
    ddg, error = get_estimate(result)
    assert pytest.approx(ddg.m) == -13.506009
    assert pytest.approx(error.m) == 0.0


def test_get_transformation_name(bace_network, bace_cleaned_result):
    cleaned_result = load_results_file(bace_cleaned_result)
    network = parse_alchemical_network(bace_network)
    transform_name = get_transform_name(cleaned_result, network)
    assert transform_name == ("complex", "spiro10", "spiro6")


@pytest.mark.parametrize("results_data, expected_output", [
    pytest.param({
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3],
        ("solvent", "lig_CHEMBL3402754_40_14", "lig_CHEMBL3402761_1_21"): [1, 2, 3]
    }, False, id="Not connected"),
    pytest.param({
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3],
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402754_40_14"): [1, 2, 3],
        ("complex", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3],
        ("complex", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402754_40_14"): [1, 2, 3]
    }, True, id="Connected"),
    pytest.param({
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2],
        ("complex", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3]
    }, False, id="Missing repeats not connected"),
    pytest.param({
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3],
        ("complex", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"): [1, 2, 3],
        ("solvent", "lig_CHEMBL3402745_200_5", "lig_CHEMBL3402754_40_14"): [1, 2, 3]
    }, True, id="Connected partial results.")
])
def test_check_is_connected(cmet_network, results_data, expected_output):
    network = parse_alchemical_network(cmet_network)
    assert check_network_is_connected(results_data=results_data, alchemical_network=network) is expected_output


# Full CLI tests
class TestScript:
    def test_cli_clean_up_not_ran(self, bace_complete_results, bace_network):
        """Make sure an error is raised if clean up was not run."""

        runner = CliRunner()
        with runner.isolated_filesystem():
            with pytest.raises(ValueError, match="has not been cleaned up please make sure you run the `results_cleanup.py` script first"):
                _ = runner.invoke(
                    gather_data,
                    [
                        "--input_alchemical_network",
                        bace_network,
                        "--output_dir",
                        "testing",
                        "--results-folder",
                        str(bace_complete_results / "results_0")
                    ],
                    catch_exceptions=False
                )

    def test_cli_missing_results(self, cmet_network, cmet_results):
        """Make sure an error is raised if the cleaned results can not be found"""
        runner = CliRunner()
        with runner.isolated_filesystem():
            with pytest.raises(FileNotFoundError, match="Can't find results directory:"):
                _ = runner.invoke(
                    gather_data,
                    [
                        "--input_alchemical_network",
                        cmet_network,
                        "--output_dir",
                        "testing",
                        "--results-folder",
                        str(cmet_results / "results_0_remove_edge")
                    ],
                    catch_exceptions=False
                )

    def test_graph_not_connected(self, bace_full_results, tmpdir):
        """Make sure an error is raised if we gather the results and have a disconnected transformation graph"""
        runner = CliRunner()
        with pytest.raises(ValueError, match="The network built from the complete results is disconnected"):
            _ = runner.invoke(
                gather_data,
                [
                    "--input_alchemical_network",
                    str(bace_full_results / "alchemicalNetwork" / "alchemical_network.json"),
                    "--output_dir",
                    str(tmpdir / "connection_test"),
                    "--results-folder",
                    str(bace_full_results / "results_0")
                ],
                catch_exceptions=False
            )

    def test_missing_fixed_network(self, bace_full_results, tmpdir):
        """Make sure an error is raised if we have extra results which are not expected"""
        runner = CliRunner()
        with pytest.raises(ValueError, match="contains a transformation which was not expected for the alchemical network"):
            _ = runner.invoke(
                gather_data,
                [
                    "--input_alchemical_network",
                    str(bace_full_results / "alchemicalNetwork" / "alchemical_network.json"),
                    "--output_dir",
                    str(tmpdir / "full_test"),
                    "--results-folder",
                    str(bace_full_results / "results_0"),
                    "--results-folder",
                    str(bace_full_results / "results_1"),
                    "--results-folder",
                    str(bace_full_results / "results_2"),
                    "--results-folder",
                    str(bace_full_results / "fixed_network" / "results_0"),
                    "--results-folder",
                    str(bace_full_results / "fixed_network" / "results_1"),
                    "--results-folder",
                    str(bace_full_results / "fixed_network" / "results_2")
                ], catch_exceptions=False
            )


    def test_full_run(self, bace_full_results, tmpdir):
        """
        Check a full successful run of the CLI and inspect the returned results.
        """
        runner = CliRunner()
        result = runner.invoke(
            gather_data,
            [
                "--input_alchemical_network",
                str(bace_full_results / "alchemicalNetwork" / "alchemical_network.json"),
                "--output_dir",
                str(tmpdir / "full_test"),
                "--results-folder",
                str(bace_full_results / "results_0"),
                "--results-folder",
                str(bace_full_results / "results_1"),
                "--results-folder",
                str(bace_full_results / "results_2")
            ]
        )
        assert result.exit_code == 0
        # make sure the missing png file warnings are printed
        assert "Can't find cleaned results file: forward_reverse_convergence.png" in result.stdout
        assert "Total results found 6/6 indicating 0 failed transformations." in result.stdout

        # make sure we have a zip folder
        archive: pathlib.Path = tmpdir / "full_test.zip"
        assert archive.exists()
        # make sure we have the folder made into the archive
        result_folder = tmpdir / "full_test"
        assert result_folder.exists()
        # load the results JSON
        result_json = json.load(
            (result_folder / "all_network_properties.json").open(mode="r"),
            cls=JSON_HANDLER.decoder
        )
        # check for the expected keys
        assert "Network_map" in result_json
        assert "transformation_scores" in result_json
        assert "ligand_scores" in result_json
        assert "DDG_estimates" in result_json

        # do some basic checks
        # check for 2 ligands and a single edge
        assert len(result_json["Network_map"]["nodes"]) == 2
        assert len(result_json["Network_map"]["edge"]) == 1
        # check the single edge has a score
        assert len(result_json["transformation_scores"]) == 1
        assert "edge_spiro2_spiro1" in result_json["transformation_scores"]
        # check the ligand scores
        assert len(result_json["ligand_scores"]) == 2
        assert "ligand_spiro2" in result_json["ligand_scores"]
        assert "ligand_spiro1" in result_json["ligand_scores"]

        # check the edge estimates
        assert len(result_json["DDG_estimates"]) == 6
        for repeat in range(3):
            for phase in ["solvent", "complex"]:
                edge_name = f"{phase}_spiro2_spiro1_repeat_{repeat}"
                assert edge_name in result_json["DDG_estimates"]
                # make sure the edge folder was created
                edge_folder = result_folder / edge_name
                assert edge_folder.exists()
                # check the files were moved to the directory
                for f_name in [
                    "structural_analysis_data.npz",
                    "energy_replica_state.npz",
                    "simulation_real_time_analysis.yaml",
                    "info.yaml"
                ]:
                    assert (edge_folder / f_name).exists()

    def test_full_run_fixed_network(self, bace_full_results, tmpdir):
        """
        Check a full successful run of the CLI with a fixed network and inspect the returned results.
        """
        runner = CliRunner()
        result = runner.invoke(
            gather_data,
            [
                "--input_alchemical_network",
                str(bace_full_results / "alchemicalNetwork" / "alchemical_network.json"),
                "--fixed_alchemical_network",
                str(bace_full_results / "fixed_network" / "alchemicalNetwork"/ "alchemical_network.json"),
                "--output_dir",
                str(tmpdir / "full_test"),
                "--results-folder",
                str(bace_full_results / "results_0"),
                "--results-folder",
                str(bace_full_results / "results_1"),
                "--results-folder",
                str(bace_full_results / "results_2"),
                "--results-folder",
                str(bace_full_results / "fixed_network" / "results_0"),
                "--results-folder",
                str(bace_full_results / "fixed_network" / "results_1"),
                "--results-folder",
                str(bace_full_results / "fixed_network" / "results_2")
            ]
        )
        assert result.exit_code == 0
        # make sure the missing png file warnings are printed
        assert "Can't find cleaned results file: forward_reverse_convergence.png" in result.stdout
        assert "Total results found 12/12 indicating 0 failed transformations." in result.stdout
        # make sure we have a zip folder
        archive: pathlib.Path = tmpdir / "full_test.zip"
        assert archive.exists()
        # make sure we have the folder made into the archive
        result_folder = tmpdir / "full_test"
        assert result_folder.exists()
        # load the results JSON
        result_json = json.load(
            (result_folder / "all_network_properties.json").open(mode="r"),
            cls=JSON_HANDLER.decoder
        )
        # check for the expected keys
        assert "Network_map" in result_json
        assert "transformation_scores" in result_json
        assert "ligand_scores" in result_json
        assert "DDG_estimates" in result_json

        # do some basic checks
        # check for 3 ligands and two edges in the combined network
        assert len(result_json["Network_map"]["nodes"]) == 3
        assert len(result_json["Network_map"]["edge"]) == 2
        # check the two edges have scores
        assert len(result_json["transformation_scores"]) == 2
        assert "edge_spiro2_spiro1" in result_json["transformation_scores"]
        assert "edge_spiro2_spiro3" in result_json["transformation_scores"]
        # check the ligand scores
        assert len(result_json["ligand_scores"]) == 3
        assert "ligand_spiro2" in result_json["ligand_scores"]
        assert "ligand_spiro1" in result_json["ligand_scores"]
        assert "ligand_spiro3" in result_json["ligand_scores"]

        # check the edge estimates
        assert len(result_json["DDG_estimates"]) == 12
        for edge in ["spiro2_spiro1", "spiro2_spiro3"]:
            for repeat in range(3):
                for phase in ["solvent", "complex"]:
                    edge_name = f"{phase}_{edge}_repeat_{repeat}"
                    assert edge_name in result_json["DDG_estimates"]
                    # make sure the edge folder was created
                    edge_folder = result_folder / edge_name
                    assert edge_folder.exists()
                    # check the files were moved to the directory
                    for f_name in [
                        "structural_analysis_data.npz",
                        "energy_replica_state.npz",
                        "simulation_real_time_analysis.yaml",
                        "info.yaml"
                    ]:
                        assert (edge_folder / f_name).exists()