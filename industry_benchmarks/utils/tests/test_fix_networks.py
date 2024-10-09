import pathlib
import pytest
from importlib import resources
import gufe
import os
import glob
import shlex

from ..fix_networks import (
    parse_alchemical_network,
    parse_results,
    alchemical_network_to_ligand_network,
    fix_network,
    cli_fix_network
)

@pytest.fixture
def results():
    with resources.files("utils.tests.data.bace_results") as d:
        yield glob.glob(f"{str(d)}/results_*_remove_edges/*json")


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


def test_parse_alchemical_network(input_alchemical_network):
    alchem_network = parse_alchemical_network(input_alchemical_network)
    assert isinstance(alchem_network, gufe.AlchemicalNetwork)
    # AlchemicalNetwork with 11 edges, meaning 22 transformations
    assert len(alchem_network.edges) == 22


def test_parse_results(results, input_alchemical_network):
    result_alchem_network = parse_results(results, input_alchemical_network)
    assert isinstance(result_alchem_network, gufe.AlchemicalNetwork)
    # only 6 edges (12 transformations) had completed successfully
    assert len(result_alchem_network.edges) == 12
    ligand_network = alchemical_network_to_ligand_network(result_alchem_network)
    assert isinstance(ligand_network, gufe.LigandNetwork)
    # The resulting LigandNetwork should have 6 edges
    assert len(ligand_network.edges) == 6


def test_fix_network(results, input_alchemical_network, output_dir, expected_transformations, tmp_path):
    temp_out_dir = tmp_path / output_dir
    fix_network(results, input_alchemical_network, temp_out_dir)
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


class TestScript:

    def test_too_few_transforms(self, cmet_network, complete_cmet_results):
        command = f"--input_alchem_network_file {cmet_network} --output_extra_transformations ./ --result_files {complete_cmet_results[0]}"
        with pytest.raises(ValueError, match="Too few transformations found for solvent_lig_CHEMBL3402745_200_5_lig_CHEMBL3402754_40_14"):
            cli_fix_network(shlex.split(command))


    def test_detect_failed_simulation(self, cmet_network, complete_cmet_results, capsys):
        command = f"--input_alchem_network_file {cmet_network} --output_extra_transformations ./ --result_files {' '.join(complete_cmet_results)}"
        with pytest.raises(ValueError, match="Too few transformations found for solvent_lig_CHEMBL3402745_200_5_lig_CHEMBL3402744_300_4"):
            cli_fix_network(shlex.split(command))

    def test_do_nothing(self):
        pass

    def test_fix_network(self, capsys, results, input_alchemical_network, output_dir, tmp_path):
        temp_out_dir = tmp_path / output_dir
        command = f"--input_alchem_network_file {input_alchemical_network} --output_extra_transformations {temp_out_dir} --result_files {' '.join(results)}"
        cli_fix_network(shlex.split(command))
        log = capsys.readouterr().out
        assert "Planned input  no. ligands: 9" in log
        assert "Planned input  no. connections: 11" in log
        assert "Simulation results no. ligands: 6" in log
        assert "Simulation results no. connections: 6" in log
        assert "Missing ligands in simulation results: 3" in log
        assert "Disconnected networks which need patching: 4" in log
        assert "Ligands in each disconnected network: [6, 1, 1, 1]" in log