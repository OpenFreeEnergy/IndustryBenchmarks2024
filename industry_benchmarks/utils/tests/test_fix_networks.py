import json
import pathlib
import pytest
from importlib import resources
import gufe
import os
import glob
import shlex

from gufe import LigandNetwork, AlchemicalNetwork
from gufe.tokenization import JSON_HANDLER

from ..fix_networks import (
    parse_alchemical_network,
    parse_results,
    alchemical_network_to_ligand_network,
    cli_fix_network,
    get_settings,
    get_settings_charge_changes,
    decompose_disconnected_ligand_network,
    get_alchemical_charge_difference,
    get_transformation_alternate
)

@pytest.fixture
def results():
    with resources.files("utils.tests.data.bace_results") as d:
        yield glob.glob(f"{str(d)}/results_*_remove_edges/*json")

@pytest.fixture
def bace_complete_results():
    with resources.files("utils.tests.data.bace_results") as d:
        yield glob.glob(f"{str(d)}/results_[0-9]/*json")


@pytest.fixture
def input_alchemical_network():
    with resources.files("utils.tests.data.bace_results") as d:
        yield str(d / "alchemicalNetwork/alchemical_network.json")


@pytest.fixture
def output_dir():
    return pathlib.Path("utils/tests/data/bace_results/newNetwork")


@pytest.fixture
def expected_transformations():
    return ["complex_spiro10_spiro1.json",
            "complex_spiro1_spiro15.json",
            "complex_spiro2_spiro15.json",
            "complex_spiro5_spiro10.json",
            "complex_spiro5_spiro6.json",
            "complex_spiro6_spiro4.json",
            "solvent_spiro10_spiro1.json",
            "solvent_spiro1_spiro15.json",
            "solvent_spiro2_spiro15.json",
            "solvent_spiro5_spiro10.json",
            "solvent_spiro5_spiro6.json",
            "solvent_spiro6_spiro4.json",
            "complex_spiro1_spiro10.json",
            "complex_spiro15_spiro1.json",
            "complex_spiro15_spiro2.json",
            "complex_spiro10_spiro5.json",
            "complex_spiro6_spiro5.json",
            "complex_spiro4_spiro6.json",
            "solvent_spiro1_spiro10.json",
            "solvent_spiro15_spiro1.json",
            "solvent_spiro15_spiro2.json",
            "solvent_spiro10_spiro5.json",
            "solvent_spiro6_spiro5.json",
            "solvent_spiro4_spiro6.json",
            "complex_spiro4_spiro10.json",
            "complex_spiro10_spiro4.json",
            "solvent_spiro4_spiro10.json",
            "solvent_spiro10_spiro4.json",
            ]


@pytest.fixture
def complete_cmet_results():
    with resources.files("utils.tests.data") as d:
        yield glob.glob(f"{str(d)}/cmet_results/results_[0-9]/*json")


@pytest.fixture
def cmet_network():
    with resources.files("utils.tests.data") as d:
        yield str(d / "cmet_results/alchemicalNetwork/alchemical_network.json")


@pytest.fixture
def eg5_network():
    with resources.files("utils.tests.data") as d:
        yield str(d / "eg5_results/subset_alchemical_network.json")

@pytest.fixture
def eg5_results():
    with resources.files("utils.tests.data") as d:
        yield glob.glob(f"{str(d)}/eg5_results/results_[0-9]/*.json")

@pytest.fixture
def bace_cleaned_result():
    with resources.files("utils.tests.data.bace_results") as d:
        yield str(d / "cleaned_results/complex_spiro10_spiro6.json")



def test_parse_alchemical_network(input_alchemical_network):
    alchem_network = parse_alchemical_network(input_alchemical_network)
    assert isinstance(alchem_network, gufe.AlchemicalNetwork)
    # AlchemicalNetwork with 11 edges, meaning 22 transformations
    assert len(alchem_network.edges) == 22


def test_parse_results_missing(results, input_alchemical_network):
    """Test parsing the results files with missing results"""
    result_alchem_network = parse_results(results, input_alchemical_network)
    assert isinstance(result_alchem_network, gufe.AlchemicalNetwork)
    # only 6 edges (12 transformations) had completed successfully
    assert len(result_alchem_network.edges) == 12
    ligand_network = alchemical_network_to_ligand_network(result_alchem_network)
    assert isinstance(ligand_network, gufe.LigandNetwork)
    # The resulting LigandNetwork should have 6 edges
    assert len(ligand_network.edges) == 6

def test_parse_results_complete(bace_complete_results, input_alchemical_network):
    """Make sure all results can be extracted for a complete network"""
    result_alchem_network = parse_results(bace_complete_results, input_alchemical_network)
    # all edges should have finished for this network
    # 11 edges with two legs should be 22 transformations
    assert len(result_alchem_network.edges) == 22
    ligand_network = alchemical_network_to_ligand_network(result_alchem_network)
    assert len(ligand_network.edges) == 11

def test_decompose_network(results, input_alchemical_network):
    """Make sure we can correctly identify the subnetworks in a failed network"""
    result_alchem_network = parse_results(results, input_alchemical_network)
    result_ligand_network = alchemical_network_to_ligand_network(result_alchem_network)
    ligand_sub_networks = decompose_disconnected_ligand_network(result_ligand_network)
    assert [len(sb.nodes) for sb in ligand_sub_networks] == [6]

def test_get_transform_alternate(input_alchemical_network, bace_cleaned_result):
    """Test we can rebuild the transform if we cleaned it up too much"""
    result = json.load(
        open(bace_cleaned_result, "r"),
        cls=JSON_HANDLER.decoder
    )
    alchemical_network = AlchemicalNetwork.from_dict(
        json.load(open(input_alchemical_network, "r"),
        cls=JSON_HANDLER.decoder
    ))
    transform, phase = get_transformation_alternate(pur=result, alchemical_network=alchemical_network)
    assert phase == "complex"



class TestScript:
    def test_too_few_transforms(self, cmet_network, complete_cmet_results):
        """Make sure an error is raised if we have too few transforms for an edge."""
        command = f"--input_alchem_network_file {cmet_network} --output_extra_transformations ./ --result_files {complete_cmet_results[0]}"
        with pytest.raises(ValueError, match="Too few transformations found for solvent_lig_CHEMBL3402745_200_5_lig_CHEMBL3402754_40_14"):
            cli_fix_network(shlex.split(command))


    def test_detect_failed_simulation(self, cmet_network, complete_cmet_results, capsys):
        """Make sure a message is printed when a simulation fails."""
        command = f"--input_alchem_network_file {cmet_network} --output_extra_transformations ./ --result_files {' '.join(complete_cmet_results)}"
        with pytest.raises(ValueError, match="Too few transformations found for solvent_lig_CHEMBL3402745_200_5_lig_CHEMBL3402744_300_4"):
            cli_fix_network(shlex.split(command))
        assert "lig_CHEMBL3402745_200_5_solvent_lig_CHEMBL3402744_300_4_solvent_solvent.json is a failed simulation" in capsys.readouterr().out

    def test_do_nothing(self, bace_complete_results, input_alchemical_network, output_dir, tmp_path, capsys):
        """Make sure we can detect when there is nothing to do to the network"""
        temp_out_dir = tmp_path / output_dir
        command = f"--input_alchem_network_file {input_alchemical_network} --output_extra_transformations {temp_out_dir} --result_files {' '.join(bace_complete_results)}"
        cli_fix_network(shlex.split(command))
        log = capsys.readouterr().out
        assert "Did not find disconnected components in alchemical network, nothing to do here!" in log

    def test_fix_network_default(self, capsys, results, input_alchemical_network, output_dir, tmp_path, expected_transformations):
        temp_out_dir = tmp_path / output_dir
        command = f"--input_alchem_network_file {input_alchemical_network} --output_extra_transformations {temp_out_dir} --result_files {' '.join(results)}"
        cli_fix_network(shlex.split(command))
        log = capsys.readouterr().out
        # make sure all expected prints are emitted
        assert "Planned input  no. ligands: 9" in log
        assert "Planned input  no. connections: 11" in log
        assert "Simulation results no. ligands: 6" in log
        assert "Simulation results no. connections: 6" in log
        assert "Missing ligands in simulation results: 3" in log
        assert "Disconnected networks which need patching: 4" in log
        assert "Ligands in each disconnected network: [6, 1, 1, 1]" in log
        assert os.path.isdir(temp_out_dir)
        output_alchemical_network = temp_out_dir / 'alchemical_network.json'
        assert output_alchemical_network.is_file()
        new_transformations = glob.glob(f"{temp_out_dir}/transformations/*json")
        # The direction of the transformation can switch, therefore include both
        # directions in the expected_transformations
        # Some of the scores of potential transformations are the same,
        # therefore the new_transformations will not always be the same
        assert len(new_transformations) == 12
        for transform in new_transformations:
            assert transform.split("/")[-1] in expected_transformations

        # load the network and check the protocol settings are right for no charge changes
        new_network = AlchemicalNetwork.from_dict(
            json.load(
                open(output_alchemical_network),
                cls=JSON_HANDLER.decoder
            )
        )
        default_settings = get_settings()
        for edge in new_network.edges:
            assert edge.protocol.settings == default_settings



    def test_fix_network_cofactor_charges(self, eg5_network, eg5_results, tmp_path, capsys):
        """Make sure we can fix networks with cofactors and charge changes."""

        temp_out_dir = tmp_path / "eg5_tranformations"
        command = f"--input_alchem_network_file {eg5_network} --output_extra_transformations {temp_out_dir} --result_files {' '.join(eg5_results)}"
        # TODO change the results to trigger the charge warning
        cli_fix_network(shlex.split(command))
        log = capsys.readouterr().out
        # make sure the network is fixed and the message is emitted
        assert "Planned input  no. ligands: 4" in log
        assert "Disconnected networks which need patching: 3" in log
        assert "Ligands in each disconnected network: [2, 1, 1]" in log
        assert "LOG: finding additional connections to merge the broken networks:" in log
        # load the networks and check the protocol settings for each leg
        output_alchemical_network = temp_out_dir / 'alchemical_network.json'
        new_network = AlchemicalNetwork.from_dict(
            json.load(
                open(output_alchemical_network),
                cls=JSON_HANDLER.decoder
            )
        )
        default_charge_settings = get_settings_charge_changes()
        defaullt_settings = get_settings()
        # make sure that charge change edges have the correct settings
        for edge in new_network.edges:
            charge_diff = get_alchemical_charge_difference(edge.mapping)
            if charge_diff != 0:
                assert edge.protocol.settings == default_charge_settings
            else:
                assert edge.protocol.settings == defaullt_settings
