import pathlib
import pytest
from importlib import resources
import gufe
import os
import glob
from ..fix_networks import (
    parse_alchemical_network,
    parse_results,
    alchemical_network_to_ligand_network,
    fix_network,
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
