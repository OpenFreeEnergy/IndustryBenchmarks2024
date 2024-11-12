from click.testing import CliRunner
import pandas as pd
import json
import pytest

from ..rename_exp_data import main as rename

@pytest.fixture()
def ligand_name_mapping() -> dict[str, str]:
    return {
        "private_ligand_1": "ligand0",
    }


@pytest.fixture(scope="function")
def csv_data() -> list[dict]:
    """Example CSV data which can be extended"""
    return [
        {
            "Ligand Name": "private_ligand_1",
            "Affinity (nM)": "105.2",
            "Affinity Error (nM)": "1.2",
            "Annotation": ""
        }
    ]

def test_full_run(tmpdir, csv_data, ligand_name_mapping):
    """
    Test a normal running of the script.
    """
    runner = CliRunner()
    pd.DataFrame(csv_data).to_csv(tmpdir / "exp_data.csv")
    with open(tmpdir / "ligand_name_mapping_PRIVATE.json", "w") as out:
        json.dump(ligand_name_mapping, out)

    result = runner.invoke(
        rename,
        [
            "--experimental-data",
            str(tmpdir / "exp_data.csv"),
            "--name-mapping-file",
            str(tmpdir / "ligand_name_mapping_PRIVATE.json"),
            "--output",
            str(tmpdir / "new_csv.csv")
        ],
    )
    assert result.exit_code == 0
    # load the new csv and check it has been renamed
    new_csv = pd.read_csv(tmpdir / "new_csv.csv")
    assert new_csv.iloc[0]["Ligand Name"] == "ligand0"


def test_missing_name_in_mapping(tmpdir, csv_data, ligand_name_mapping):
    """
    Make sure an error is raised if we can not convert a name in the exp csv file.
    """

    runner = CliRunner()

    # add an entry which is not mapped
    csv_data.append({
        "Ligand Name": "missing_ligand_1",
        "Affinity (nM)": "200.2",
        "Affinity Error (nM)": "2.2",
        "Annotation": "Racemate"
    })
    pd.DataFrame(csv_data).to_csv(tmpdir / "exp_data.csv")

    with open(tmpdir / "ligand_name_mapping_PRIVATE.json", "w") as out:
        json.dump(ligand_name_mapping, out)


    with pytest.raises(RuntimeError, match="Could not convert missing_ligand_1 as it was not found in the name mapping."):
        _ = runner.invoke(
            rename,
            [
                "--experimental-data",
                str(tmpdir / "exp_data.csv"),
                "--name-mapping-file",
                str(tmpdir / "ligand_name_mapping_PRIVATE.json"),
                "--output",
                str(tmpdir / "new_csv.csv")
            ],
            catch_exceptions=False
        )
