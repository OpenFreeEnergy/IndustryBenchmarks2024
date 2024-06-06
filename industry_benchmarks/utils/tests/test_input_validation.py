import pytest
from importlib import resources
import click.testing
from ..input_validation import run_inputs as main


@pytest.fixture
def protein():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield str(d / 'protein.pdb')


@pytest.fixture
def bad_protein():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield str(d / 'bad_protein.pdb')


@pytest.fixture
def protein_with_ptm():
    with resources.files('utils.tests.data.thrombin_inputs') as d:
        yield str(d / 'protein_TYS.pdb')


@pytest.fixture
def protein_convert_ptm():
    with resources.files('utils.tests.data.thrombin_inputs') as d:
        yield str(d / 'protein.pdb')


@pytest.fixture
def cofactors():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield str(d / 'cofactors.sdf')


class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: run-inputs" in result.output

    def test_protein(self, protein):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ['--pdb', protein])
            assert result.exit_code == 0
            assert "SIMULATION COMPLETE" in result.output

    def test_protein_cofactors(self, protein, cofactors):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                main,
                ['--pdb', protein, '--cofactors', cofactors]
            )
            assert result.exit_code == 0
            assert "SIMULATION COMPLETE" in result.output

    def test_bad_protein(self, bad_protein):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ['--pdb', bad_protein])
            assert result.exit_code == 1
            assert "SIMULATION COMPLETE" not in result.output

    def test_ptm_protein_fail(self, protein_with_ptm):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ['--pdb', protein_with_ptm])
            print(result)
            assert result.exit_code == 1
            assert "SIMULATION COMPLETE" not in result.output

    def test_revert_ptm_protein_pass(self, protein_convert_ptm):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ['--pdb', protein_convert_ptm])
            print(result)
            assert result.exit_code == 0
            assert "SIMULATION COMPLETE" in result.output
