import pytest
from importlib import resources
import os
import csv
import click.testing
from ..extras.extract_results import extract as main


@pytest.fixture
def results_0():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield d / 'results_0'

@pytest.fixture
def results_1():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield d / 'results_1'

@pytest.fixture
def results_2():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield d / 'results_2'

@pytest.fixture
def outfile():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield d / 'ddg.tsv'

@pytest.fixture
def outfile_dG():
    with resources.files('utils.tests.data.bace_inputs') as d:
        yield d / 'dg.tsv'

class TestScript:
    def test_invoke(self):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(main, ["--help"])
            assert result.exit_code == 0
            assert "Usage: extract" in result.output

    def test_extract_results(self, results_0, results_1, results_2, outfile, outfile_dG):
        runner = click.testing.CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                main,
                ['--results_0', results_0, '--results_1', results_1,
                 '--results_2', results_2, '-o', outfile,
                 '-o_DG', outfile_dG]
            )
            assert result.exit_code == 0
            assert os.path.exists(str(outfile))
            assert os.path.exists(str(outfile_dG))

            # Check the DDG file
            with open(str(outfile), 'r') as f:
                ddg = csv.reader(f, delimiter=',', quotechar='"')
                headers = next(ddg)
                edges = [r for r in ddg]
                # Results from 11 edges are present in the file
                assert len(edges) == 11

            # Check the DG file
            with open(str(outfile_dG), 'r') as f:
                dg = csv.reader(f, delimiter=',', quotechar='"')
                headers = next(dg)
                nodes = [r for r in dg]
                # Results from 9 ligands are present in the file
                assert len(nodes) == 9

