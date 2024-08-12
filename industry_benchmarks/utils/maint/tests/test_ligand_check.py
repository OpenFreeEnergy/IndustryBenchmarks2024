import pytest
from importlib import resources
import click.testing
from ..check_ligand_prep import run as main


@pytest.fixture
def dlk_ligands():
    with resources.files('utils.maint.data') as d:
        yield str(d / 'dlk_ligands.sdf')


@pytest.fixture
def hif2a_ligands():
    with resources.files('utils.maint.data') as d:
        yield str(d / 'hif2a_automap_ligands.sdf')


@pytest.fixture
def tnks2_ligands():
    with resources.files('utils.maint.data') as d:
        yield str(d / 'tnks2_fullmap_ligands.sdf')


@pytest.fixture
def tyk2_ligands():
    with resources.files('utils.maint.data') as d:
        yield str(d / 'tyk2_ligands.sdf')


class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: run" in result.output

    @pytest.mark.parametrize('ligands,num', [
        ["dlk_ligands", 2],
        ["hif2a_ligands", 1],
        ["tnks2_ligands", 6],
        ["tyk2_ligands", 2],
    ])
    def test_ligands(self, request, ligands, num):
        ligands = request.getfixturevalue(ligands)
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ['--ligands', ligands])
            assert result.exit_code == 0
            assert f"Error: Number of duplicate ligands found:  {num}" in result.output
