import pathlib

import pytest
from importlib import resources
import os
import json
from gufe.tokenization import JSON_HANDLER
import click.testing
from openfe.setup import LigandNetwork
from openfe.protocols import openmm_rfe
from ..plan_rbfe_network import run_inputs as main
from ..plan_rbfe_network import get_settings_charge_changes


@pytest.fixture
def ligands():
    with resources.files('utils.tests.data.eg5_inputs') as d:
        yield str(d / 'ligands_subset.sdf')


@pytest.fixture
def protein():
    with resources.files('utils.tests.data.eg5_inputs') as d:
        yield str(d / 'protein.pdb')


@pytest.fixture
def cofactors():
    with resources.files('utils.tests.data.eg5_inputs') as d:
        yield str(d / 'cofactor.sdf')


@pytest.fixture
def output():
    return pathlib.Path('utils.tests.data.eg5_inputs.input_jsons')


class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: run-inputs" in result.output

    def test_run_inputs(self, ligands, protein, cofactors, output):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            # Test warning for charge changing transformations
            with pytest.warns(UserWarning, match='Charge changing transformation'):
                result = runner.invoke(
                    main,
                    ['--ligands', ligands, '--pdb', protein,
                     '--cofactors', cofactors, '--output', output]
                )
 
                with open(output / "ligand_network.graphml") as f:
                    graphml = f.read()
                ligand_network = LigandNetwork.from_graphml(graphml)
                for edge in ligand_network.edges:
                    assert edge.componentA.to_openff().partial_charges is not None
                output_files = list(output.glob("*.json"))
                assert result.exit_code == 0
                assert os.path.exists(str(output))
                # Check if the output files are created
                assert len(output_files) == 6
 
                # Check that we have 4 json files with charge change settings, 
                # 2 json files with default settings
                charge_change_json = []
                default_json = []
                for f in output_files:
                    d = json.load(open(f, 'r'), cls=JSON_HANDLER.decoder)
                    protocol = openmm_rfe.RelativeHybridTopologyProtocol.from_dict(d['protocol'])
                    lambda_windows = protocol.settings.lambda_settings.lambda_windows
                    if lambda_windows == 22:
                        charge_change_json.append(f)
                    if lambda_windows == 11:
                        default_json.append(f)
                assert len(charge_change_json) == 4
                assert len(default_json) == 2
