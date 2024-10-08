import pytest
from importlib import resources
import os
import csv
import click.testing
from ..extras.extract_results import extract as main


@pytest.fixture
def results_0():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_0"


@pytest.fixture
def results_0_partial():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_0_partial"


@pytest.fixture
def results_0_remove_edge():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_0_remove_edge"


@pytest.fixture
def results_1():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_1"


@pytest.fixture
def results_1_remove_edge():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_1_remove_edge"


@pytest.fixture
def results_2():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_2"


@pytest.fixture
def results_2_remove_edge():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/results_2_remove_edge"


@pytest.fixture
def ligand_network():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/alchemicalNetwork/ligand_network.graphml"


@pytest.fixture
def outfile():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/ddg.tsv"


@pytest.fixture
def outfile_dG():
    with resources.files("utils.tests.data") as d:
        yield d / "cmet_results/dg.tsv"


class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: extract" in result.output

    def test_extract_results(
        self, results_0, results_1, results_2, ligand_network, outfile, outfile_dG
    ):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                main,
                [
                    "--results_0",
                    results_0,
                    "--results_1",
                    results_1,
                    "--results_2",
                    results_2,
                    "--input_ligand_network_file",
                    ligand_network,
                    "-o",
                    outfile,
                    "-o_DG",
                    outfile_dG,
                ],
            )
            assert result.exit_code == 0
            assert os.path.exists(str(outfile))
            assert os.path.exists(str(outfile_dG))

            # Check the DDG file
            with open(str(outfile), "r") as f:
                ddg = csv.reader(f, delimiter=",", quotechar='"')
                headers = next(ddg)
                edges = [r for r in ddg]
                # Results from 4 edges are present in the file
                assert len(edges) == 4

            # Check the DG file
            with open(str(outfile_dG), "r") as f:
                dg = csv.reader(f, delimiter=",", quotechar='"')
                headers = next(dg)
                nodes = [r for r in dg]
                # Results from 5 ligands are present in the file
                assert len(nodes) == 5


    def test_extract_results_removed_edge(
        self,
        results_0_remove_edge,
        results_1_remove_edge,
        results_2_remove_edge,
        ligand_network,
        outfile,
        outfile_dG
    ):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                main,
                [
                    "--results_0",
                    results_0_remove_edge,
                    "--results_1",
                    results_1_remove_edge,
                    "--results_2",
                    results_2_remove_edge,
                    "--input_ligand_network_file",
                    ligand_network,
                    "-o",
                    outfile,
                    "-o_DG",
                    outfile_dG,
                ],
            )
            assert result.exit_code == 0
            assert os.path.exists(str(outfile))
            assert os.path.exists(str(outfile_dG))

            # Check the DDG file
            with open(str(outfile), "r") as f:
                ddg = csv.reader(f, delimiter=",", quotechar='"')
                headers = next(ddg)
                edges = [r for r in ddg]
                # Results from 3 edges are present in the file
                assert len(edges) == 3

            # Check the DG file
            with open(str(outfile_dG), "r") as f:
                dg = csv.reader(f, delimiter=",", quotechar='"')
                headers = next(dg)
                nodes = [r for r in dg]
                # Results from 4 ligands are present in the file
                assert len(nodes) == 4


    def test_extract_results_partial(
        self, results_0_partial, results_1, results_2, ligand_network, outfile, outfile_dG
    ):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                main,
                [
                    "--results_0",
                    results_0_partial,
                    "--results_1",
                    results_1,
                    "--results_2",
                    results_2,
                    "--input_ligand_network_file",
                    ligand_network,
                    "-o",
                    outfile,
                    "-o_DG",
                    outfile_dG,
                ],
            )
            assert result.exit_code == 1
            assert "Some files are partially missing for some edges" in result.exception.args[0]
