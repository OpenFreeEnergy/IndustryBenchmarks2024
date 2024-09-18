import pathlib

import pytest
from importlib import resources
import os
import json
import glob
from gufe.tokenization import JSON_HANDLER
import click.testing
from openfe.setup import LigandNetwork
from openfe.protocols import openmm_rfe
from ..tape_networks import cli_do_taping as main


@pytest.fixture
def results():
    with resources.files("utils.tests.data.bace_inputs") as d:
        yield pathlib.Path(d / "results_broken_network/*json")


@pytest.fixture
def input_alchemical_network():
    with resources.files("utils.tests.data.bace_inputs") as d:
        yield str(d / "alchemical_network.json")


@pytest.fixture
def output_dir():
    return pathlib.Path("utils.tests.data.bace_inputs.tapesForAlchemicalNetwork")


class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: cli-do-taping" in result.output

    def test_do_taping(self, results, input_alchemical_network, output_dir):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            # Test warning for charge changing transformations
            result = runner.invoke(
                main,
                ["-rp", results, "-ia", input_alchemical_network, "-oa", output_dir],
            )
            assert result.exit_code == 0
            assert os.path.exists(str(output_dir))
            # Create a list of all .json files
            output_files = list(output_dir.glob("transformations/*.json"))
            # Check if the output files are created
            assert len(output_files) == 12
            for f in output_files:
                assert f.is_file()
