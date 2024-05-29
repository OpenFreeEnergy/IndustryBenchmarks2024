import pathlib

import pytest
from importlib import resources
import os
import click.testing
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
    return pathlib.Path('input_jsons')



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
            result = runner.invoke(
                main,
                ['--ligands', ligands, '--pdb', protein,
                 '--cofactors', cofactors, '--output', output]
            )
            output_files = list(output.glob("*.json"))
            assert result.exit_code == 0
            assert os.path.isdir(output)
            # Check if the output files are created
            assert len(output_files) > 0

    def test_get_settings_charge_changes(self):
        get_settings_charge_changes()


    # Test that the number of transformation files for solvent and complex are equal
    #


