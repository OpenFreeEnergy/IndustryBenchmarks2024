import click
from cinnabar.plotting import _master_plot, plot_DGs
from cinnabar import FEMap
import pandas as pd
import pathlib
from openff.units import unit
import json
from gufe.tokenization import JSON_HANDLER
import numpy as np
import tqdm
from typing import Literal
import pymbar
from pymbar import MBAR
import shutil
from itertools import combinations


def load_exp_data(filename: pathlib.Path) -> pd.DataFrame:
    """
    Load the experimental dataframe and perform some checks on the data.

    Parameters
    ----------
    filename: pathlib.Path
        The path to the experimental data file.

    Returns
    -------
        The extracted dataframe which has been validated.

    Raises
    ------
    ValueError: if a ligand name appears in the csv file more than once or the column names are not as expected.

    Notes
    -----
    Any missing error values are changed to 0.0
    """
    expected_names = [
        "Ligand Name",
        "Affinity (nM)",
        "Affinity Error (nM)",
        "Annotation",
    ]
    exp_data = pd.read_csv(filename.as_posix())
    # check the column names match what we expect
    for column in exp_data.columns:
        if column not in expected_names:
            raise ValueError(
                f"Column {column} was not expected check the naming matches the specification."
            )
    unique_names = exp_data["Ligand Name"].unique()
    if len(unique_names) != len(exp_data):
        raise ValueError(
            "There are repeated ligand names in the experimental data CSV file."
        )
    # do something to clean non-float values in the affinity or error columns
    return exp_data

def get_exp_ddg_fep_plus(
    experimental_data: pd.DataFrame, ligand_a: str, ligand_b: str
) -> tuple[float, float]:
    """
    Given the experimental dataframe reported by fep plus and the ligand a & b identifiers get the reference DDG and
    error.

    Parameters
    ----------
    experimental_data: pd.DataFrame
        The experimental dataframe reported by FEP plus
    ligand_a, ligand_b: str
        The identifiers of the ligands involved in the edge, computed as A to B.
    Returns
    -------
        The experimental DDG and error estimate for this transformation.
    """
    ligand_a_dg = experimental_data[
        experimental_data["Ligand name"] == ligand_a
        ].iloc[0]["Exp. dG (kcal/mol)"]
    ligand_b_dg = experimental_data[
        experimental_data["Ligand name"] == ligand_b
        ].iloc[0]["Exp. dG (kcal/mol)"]
    ddg = ligand_b_dg - ligand_a_dg
    return ddg, 0.0


def get_exp_ddg(
    experimental_data: pd.DataFrame, ligand_a: str, ligand_b: str
) -> tuple[float, float]:
    """
    Given the experimental dataframe and the ligand a & b identifiers get the reference DDG and error.

    Assumes the transform is from A to B and returns DDG as B-A.

    Parameters
    ----------
    experimental_data: pd.DataFrame
        The experimental dataframe in the format specified in the industry benchmarking docs.
    ligand_a, ligand_b: str
        The identifiers of the ligands involved in the edge, computed as A to B.

    Returns
    -------
        The experimental DDG and error estimate for this transformation.

    Notes
    -----
    We assume all error validation has already been completed before running this function.
    """
    ligand_a_data = experimental_data[
        experimental_data["Ligand Name"] == ligand_a
    ].iloc[0]
    ligand_b_data = experimental_data[
        experimental_data["Ligand Name"] == ligand_b
    ].iloc[0]
    # convert from ki to dg
    ligand_a_dg, ligand_a_dg_error = ki_to_dg(
        ki=ligand_a_data["Affinity (nM)"] * unit.nanomolar,
        uncertainty=ligand_a_data["Affinity Error (nM)"] * unit.nanomolar,
    )
    ligand_b_dg, ligand_b_dg_error = ki_to_dg(
        ki=ligand_b_data["Affinity (nM)"] * unit.nanomolar,
        uncertainty=ligand_b_data["Affinity Error (nM)"] * unit.nanomolar,
    )
    ddg = ligand_b_dg - ligand_a_dg
    delta_ddg = np.sqrt(ligand_a_dg_error**2 + ligand_b_dg_error**2)
    return ddg.m, delta_ddg.m


def ki_to_dg(
    ki: unit.Quantity,
    uncertainty: unit.Quantity,
    temperature: unit.Quantity = 298.15 * unit.kelvin,
) -> tuple[unit.Quantity, unit.Quantity]:
    """
    Convenience method to convert a Ki w/ a given uncertainty to an
    experimental estimate of the binding free energy.

    Args:
        ki: Experimental Ki value (e.g. 5 * unit.nanomolar)
        uncertainty: Experimental error. Note: returns 0 if =< 0 * unit.nanomolar.
        temperature: Experimental temperature. Default: 298.15 * unit.kelvin.

    Returns:
        dg: Gibbs binding free energy.
        ddg: Error in binding free energy.
    """
    if ki > 1e-15 * unit.molar:
        dg = (
            unit.molar_gas_constant
            * temperature.to(unit.kelvin)
            * np.log(ki / unit.molar)
        ).to(unit.kilocalorie_per_mole)
    else:
        raise ValueError("negative Ki values are not supported")
    # propagate the uncertainty <https://en.wikipedia.org/wiki/Propagation_of_uncertainty>
    if uncertainty > 0 * unit.molar:
        ddg = (
            unit.molar_gas_constant * temperature.to(unit.kelvin) * uncertainty / ki
        ).to(unit.kilocalorie_per_mole)
    else:
        ddg = 0 * unit.kilocalorie_per_mole

    return dg, ddg


def extract_general_edge_data(
    all_properties: dict[str, dict],
    dataset_name: str,
    industry_partner: str,
    experimental_data: pd.DataFrame,
    archive: pathlib.Path,
    is_public: bool
) -> pd.DataFrame:
    """
    Extract the top level edge information from the all properties JSON file made by the data gathering script.

    Parameters
    ----------
    all_properties: dict
        A dictionary of all network properties created by the data_gathering.py script.
    dataset_name: str
        The unique name which should be used to identify this dataset.
    industry_partner: str
        The unique name which should be used to identify the partner this dataset came from.
    experimental_data: pd.DataFrame
        The experimental CSV file from which we should extract the DDG values.
    archive: pathlib.Path
        The location of the archive to extract edge matrix data from.
    is_public: bool
        If the dataset is from the public benchmark, this changes how the exp data is extracted.

    Returns
    -------
        A DataFrame of the edge data, such as edge scores and DG values per phase and repeat.
    """
    all_data = []
    for edge_data in tqdm.tqdm(
        all_properties["Edges"].values(),
        desc=f"Extracting edge data for {dataset_name}",
        total=len(all_properties["Edges"]),
    ):
        # edge data has all the edge scores already
        # and some identifiers to group the edges later
        edge_data["dataset_name"] = dataset_name
        edge_data["partner_id"] = industry_partner
        # add in the experimental estimate for the edge
        if is_public:
            ddg, error = get_exp_ddg_fep_plus(
                experimental_data=experimental_data,
                ligand_a=edge_data["ligand_A"],
                ligand_b=edge_data["ligand_B"]
            )
        else:
            ddg, error = get_exp_ddg(
                experimental_data=experimental_data,
                ligand_a=edge_data["ligand_A"],
                ligand_b=edge_data["ligand_B"],
            )
        edge_data["exp DDG (kcal/mol)"] = ddg
        edge_data["exp dDDG (kcal/mol)"] = error
        # add in the per repeat value for each phase
        for phase in ["solvent", "complex"]:
            for repeat in range(3):
                repeat_name = f"{phase}_{edge_data['ligand_A']}_{edge_data['ligand_B']}_repeat_{repeat}"
                phase_repeat_name = f"{phase}_repeat_{repeat}"
                dg, error = all_properties["DDG_estimates"][repeat_name]
                edge_data[f"{phase_repeat_name}_DG (kcal/mol)"] = dg.to(
                    unit.kilocalorie_per_mole
                ).m
                edge_data[f"{phase_repeat_name}_dDG (kcal/mol)"] = error.to(
                    unit.kilocalorie_per_mole
                ).m
                # extract maximum com drift and ligand RMSD values
                com_drift_max = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="com"
                ).max()
                rmsd_max = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="ligand_rmsd"
                ).max()
                edge_data[f"{phase_repeat_name}_com_drift_max"] = com_drift_max
                edge_data[f"{phase_repeat_name}_ligand_rmsd_max"] = rmsd_max
                # extract the minimum off diagonal in the overlap matrix
                potential = get_edge_matrix(
                    edge_name=repeat_name,
                    archive=archive,
                    matrix_type="reduced_potential",
                )
                samples = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="samples"
                )
                # transpose the matrix as it was formated for readability
                overlap_matrix = compute_overlap_matrix(
                    reduced_potential=potential.T, samples=samples
                )
                smallest_off_diagonal = get_smallest_off_diagonal(overlap_matrix)
                edge_data[f"{phase_repeat_name}_smallest_overlap"] = (
                    smallest_off_diagonal
                )

        # once we have every repeat for each phase add the row
        all_data.append(edge_data)

    return pd.DataFrame(all_data)


def download_zenodo_archive(address: str, output_folder: pathlib.Path) -> pathlib.Path:
    """
    Download the zenodo archive to the local output folder.

    Parameters
    ----------
    address: str
        The http address to the zenodo archive.
    output_folder: pathlib.Path
        The path to the output folder the data should be downloaded to.

    Returns
    -------
    The path to the downloaded archive
    """
    import requests
    import wget
    zenodo_records_url = "https://zenodo.org/api/records/"
    # get the record ID
    if "records" in address:
        record_id = address.split("/")[-1]
        record_id = record_id.strip()
        # get the data for the record should be a single zip folder
        record_data = requests.get(zenodo_records_url + record_id)

        if record_data.ok:
            record_json = record_data.json()
            # for each of the files in the folder download them to the output folder
            record_files = [f for f in record_json["files"]]
            # should be a single folder
            assert len(record_files) == 1, print("More than a single folder found!")
            files = []
            for f in record_files:
                filename = f.get("filename") or f["key"]
                file_address = f"https://zenodo.org/records/{record_id}/files/{filename}"
                wget.download(file_address, output_folder.as_posix())
                files.append(filename)
            return output_folder.joinpath(files[0])
        else:
            record_data.raise_for_status()

    else:
        raise RuntimeError(f"The zenodo address ({address}) did not contain the record id")



def get_edge_matrix(
    edge_name: str,
    archive: pathlib.Path,
    matrix_type: Literal[
        "com", "ligand_rmsd", "protein_rmsd", "reduced_potential", "samples"
    ],
) -> np.array:
    """
    Given the edge name for this run and the archive extract the matrix file.

    Parameters
    ----------
    edge_name: str
        The name of the edge of the format phase_ligandA_ligand_B_repeat_x
    archive: pathlib.Path
        The path to the extracted zenodo archive.
    matrix_type: str
        The shorthand for the type of matrix to extract.

    Returns
    -------
    The loaded matrix of the format (n_samples, n_replicas)
    """
    types_to_files = {
        "com": "ligand_wander.txt",
        "ligand_rmsd": "ligand_RMSD.txt",
        "protein_rmsd": "protein_RMSD.txt",
        "reduced_potential": "u_ln.txt",
        "samples": "N_l.txt",
    }
    file_name = archive.joinpath(edge_name, types_to_files[matrix_type])
    return np.loadtxt(file_name)


def compute_overlap_matrix(reduced_potential: np.array, samples: np.array) -> np.array:
    """
    Compute the MBAR overlap matrix using pymbar3.

    Parameters
    ----------
    reduced_potential: np.array
        The reduced potential calculated for each sample in the form (n_states, n_samples)
    samples
        The number of samples per state.
    Returns
    -------
        The MBAR overlap matrix calculated with pymbar.
    """
    mbar = MBAR(reduced_potential, samples)
    _ = mbar.getFreeEnergyDifferences()
    overlap_matrix = mbar.computeOverlap()
    return overlap_matrix["matrix"]


def get_smallest_off_diagonal(matrix: np.array) -> float:
    """
    For the symmetric matrix return the smallest of diagonal value.

    Parameters
    ----------
    matrix: np.array
        The matrix we want the smallest of diagonal value for.

    Returns
    -------
    """
    diagonal = np.diagonal(matrix, offset=1)
    return diagonal.min()

def get_cumulative_dg(
        reduced_potential: np.array,
        samples: np.array
) -> dict[str, float]:
    """
    Compute the cumulative DG estimates using a % of the data from 10 to 100.

    Parameters
    ----------
    reduced_potential: np.array
        The reduced potential calculated for each sample in the form (n_states, n_samples).
    samples: np.array
        The number of samples per state.

    Returns
    -------
    A dict of the % of samples and the estimated DG and dDG in kcal/mol.
    """
    cumulative_dgs = {}
    # assuming we ran at the default openfe temperature in the protocol
    kt = unit.molar_gas_constant * 298.15 * unit.kelvin
    for i in range(1, 11):
        sliced_samples = np.floor(samples * (i / 10))
        mbar = MBAR(reduced_potential[:,:int(sliced_samples[0] * 11)], sliced_samples)
        try:
            estimate = mbar.getFreeEnergyDifferences(return_dict=True)
            DG, dDG = estimate["Delta_f"][0, -1] * kt, estimate["dDelta_f"][0, -1] * kt
            cumulative_dgs[f"Samples {int(i * 10)}% DG"] = DG.to(unit.kilocalorie_per_mole).m
            cumulative_dgs[f"Samples {int(i * 10)}% dDG"] = dDG.to(unit.kilocalorie_per_mole).m
        except pymbar.utils.ParameterError:
            # this is raised if we have an error solving
            # in this case just return na?
            cumulative_dgs[f"Samples {int(i * 10)}% DG"] = pd.NA
            cumulative_dgs[f"Samples {int(i * 10)}% dDG"] = pd.NA
    return cumulative_dgs

def extract_cumulative_dg_estimates(
    all_properties: dict[str, dict],
    dataset_name: str,
    industry_partner: str,
    archive: pathlib.Path
) -> pd.DataFrame:
    """
    For each edge in the dataset estimate the DG cumulative values using 10%, 20%, ... 100% of the data.

    Parameters
    ----------
     all_properties: dict
        A dictionary of all network properties created by the data_gathering.py script.
    dataset_name: str
        The unique name which should be used to identify this dataset.
    industry_partner: str
        The unique name which should be used to identify the partner this dataset came from.
    archive: pathlib.Path
        The location of the archive to extract edge matrix data from.

    Returns
    -------
        A dataframe of DG estimates using portions of the total data, where each row is a unique phase and repeat.
    """
    all_data = []
    for edge_data in tqdm.tqdm(
            all_properties["Edges"].values(),
            desc=f"Calculating cumulative edge data for {dataset_name}",
            total=len(all_properties["Edges"]),
    ):
        # add in the per repeat value for each phase
        for phase in ["solvent", "complex"]:
            for repeat in range(3):
                repeat_name = f"{phase}_{edge_data['ligand_A']}_{edge_data['ligand_B']}_repeat_{repeat}"
                # collect some general data
                dg_data = {
                    "ligand_A": edge_data["ligand_A"],
                    "ligand_B": edge_data["ligand_B"],
                    "phase": phase,
                    "repeat": repeat,
                    "partner_id": industry_partner
                }
                potential = get_edge_matrix(
                    edge_name=repeat_name,
                    archive=archive,
                    matrix_type="reduced_potential",
                )
                samples = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="samples"
                )
                cumulative_dg = get_cumulative_dg(reduced_potential=potential.T, samples=samples)
                # merge the cumulative data
                dg_data.update(cumulative_dg)

                all_data.append(dg_data)

    return pd.DataFrame(all_data)

def plot_ddg_vs_experiment(edge_dataframe: pd.DataFrame, output: pathlib.Path, dataset_name: str):
    """
    For the given dataframe of general edge data plot the DDG vs the experimental data using cinnabar.

    Parameters
    ----------
    edge_dataframe: pd.DataFrame
        The dataframe with general edge data generated using extract_general_edge_data
    output: pathlib.Path
        The path of the output directory, the file will be called DDG_vs_exp.svg
    dataset_name: str
        The name of the dataset which will be added as a plot title.
    """

    exp_ddg = edge_dataframe["exp DDG (kcal/mol)"]
    exp_error = edge_dataframe["exp dDDG (kcal/mol)"]
    # calculate the DDG using the std between repeats as the error
    complex_data = edge_dataframe[[f"complex_repeat_{i}_DG (kcal/mol)" for i in range(3)]]
    complex_dg = complex_data.mean(axis=1)
    complex_error = complex_data.std(axis=1)
    solvent_data = edge_dataframe[[f"solvent_repeat_{i}_DG (kcal/mol)" for i in range(3)]]
    solvent_dg = solvent_data.mean(axis=1)
    solvent_error = solvent_data.std(axis=1)
    ddg = complex_dg - solvent_dg
    # propagate errors
    ddg_error = (complex_error ** 2 + solvent_error ** 2) ** 0.5
    # create the plot
    _master_plot(
        x=exp_ddg,
        xerr=exp_error,
        y=ddg,
        yerr=ddg_error,
        statistic_type="mle",
        title=dataset_name,
        figsize=4,
        filename=output.joinpath("DDG_vs_exp.svg").as_posix()
    )

def plot_ddg_repeats(edge_dataframe: pd.DataFrame, output: pathlib.Path):
    """
    For the edge dataframe plot estimated DDG from each repeat against each other using cinnabar.

    Parameters
    ----------
    edge_dataframe: pd.DataFrame
        The dataframe with general edge data generated using extract_general_edge_data
    output: pathlib.Path
        The path of the output directory, the files will be saved as DDG_0_vs_1.svg, DDG_1_vs_2.svg, DDG_0_vs_2.svg
    """
    for repeat1, repeat2 in combinations(range(3), 2):
        repeat1_ddg = edge_dataframe[f"complex_repeat_{repeat1}_DG (kcal/mol)"] - edge_dataframe[f"solvent_repeat_{repeat1}_DG (kcal/mol)"]
        repeat1_error = (
            edge_dataframe[f"complex_repeat_{repeat1}_dDG (kcal/mol)"] ** 2 + edge_dataframe[f"solvent_repeat_{repeat1}_dDG (kcal/mol)"] ** 2
        ) ** 0.5
        repeat2_ddg = edge_dataframe[f"complex_repeat_{repeat2}_DG (kcal/mol)"] - edge_dataframe[f"solvent_repeat_{repeat2}_DG (kcal/mol)"]
        repeat2_error = (
            edge_dataframe[f"complex_repeat_{repeat2}_dDG (kcal/mol)"] ** 2 + edge_dataframe[f"solvent_repeat_{repeat2}_dDG (kcal/mol)"] ** 2
        ) ** 0.5
        _master_plot(
            x=repeat1_ddg,
            xerr=repeat1_error,
            y=repeat2_ddg,
            yerr=repeat2_error,
            statistic_type="mle",
            title="DDG between repeats",
            figsize=4,
            filename=output.joinpath(f"DDG_{repeat1}_vs_{repeat2}.svg").as_posix(),
            xlabel=f"Repeat {repeat1}",
            ylabel=f"Repeat {repeat2}"
        )

def calculate_and_plot_dg(
        edge_dataframe: pd.DataFrame,
        output: pathlib.Path,
        experimental_data: pd.DataFrame,
        is_public: bool,
        dataset_name: str
) -> pd.DataFrame:
    """
    Calculate the DG values based on the edge dataframe and plot them against the experimental data.

    Parameters
    ----------
    edge_dataframe: pd.DataFrame
        The edge data extracted by extract_general_edge_dat.
    output: pathlib.Path
        The path of the output directory the plot should be saved to (DG_vs_exp.svg)
    experimental_data: pd.DataFrame
        The experimental dataframe.
    is_public: bool
        If the dataset is from the public benchmark, this changes how the exp data is extracted.
    dataset_name: str
        The name of the dataset which should be used as the title.

    Returns
    -------
        A dataframe of the calculated DG values compared with experimental values.
    """
    # calculate the DDG using the std between repeats as the error
    complex_data = edge_dataframe[[f"complex_repeat_{i}_DG (kcal/mol)" for i in range(3)]]
    complex_dg = complex_data.mean(axis=1)
    complex_error = complex_data.std(axis=1)
    solvent_data = edge_dataframe[[f"solvent_repeat_{i}_DG (kcal/mol)" for i in range(3)]]
    solvent_dg = solvent_data.mean(axis=1)
    solvent_error = solvent_data.std(axis=1)
    ddg = complex_dg - solvent_dg
    # propagate errors
    ddg_error = (complex_error ** 2 + solvent_error ** 2) ** 0.5
    fe_map = FEMap()
    # add each prediction
    for index, row in edge_dataframe.iterrows():
        fe_map.add_relative_calculation(
            labelA=row["ligand_A"],
            labelB=row["ligand_B"],
            value=ddg.loc[index],
            uncertainty=ddg_error.loc[index],
            source="calculated",
        )
    # add the exp data points
    for _, row in experimental_data.iterrows():
        if is_public:
            fe_map.add_experimental_measurement(
                label=row["Ligand name"],
                value=row["Exp. dG (kcal/mol)"] * unit.kilocalorie_per_mole,
                uncertainty=0.0 * unit.kilocalorie_per_mole
            )
        else:
            fe_map.add_experimental_measurement(
                label=row["Ligand Name"],
                value=row["Affinity (nM)"] * unit.nanomolar,
                uncertainty=row["Affinity Error (nM)"] * unit.nanomolar
            )

    fe_map.generate_absolute_values()
    plot_DGs(
        graph=fe_map.to_legacy_graph(),
        title=dataset_name,
        filename=output.joinpath("DG_vs_exp.svg").as_posix(),
        figsize=4
    )
    # if we have a public result plot DG vs FEP+
    absolute_df = fe_map.get_absolute_dataframe()
    abs_calc = absolute_df[absolute_df["computational"] == True]
    if is_public:
        openfe, openfe_error, fep_plus, fep_plus_error = [], [], [], []
        for _, row in abs_calc.iterrows():
            openfe.append(row["DG (kcal/mol)"])
            openfe_error.append(row["uncertainty (kcal/mol)"])
            fep_plus_data = experimental_data[experimental_data["Ligand name"] == row["label"]].iloc[0]
            fep_plus.append(fep_plus_data["Pred. dG (kcal/mol)"])
            fep_plus_error.append(fep_plus_data["Pred. dG std. error (kcal/mol)"])
        # plot the data
        _master_plot(
            # shift by the mean Exp value to match the FEP+ data
            x=np.array(openfe) + experimental_data["Exp. dG (kcal/mol)"].mean(axis=0),
            xerr=np.array(openfe_error),
            y=np.array(fep_plus),
            yerr=np.array(fep_plus_error),
            title="OpenFE vs FEP+",
            xlabel="OpenFE",
            ylabel="FEP+",
            quantity=rf"$\Delta$ G",
            filename=output.joinpath("openfe_vs_fep+.svg").as_posix()
        )

    return abs_calc.drop(columns=["computational", "source"])



@click.command()
@click.option(
    "-z",
    "--zenodo",
    help="The address of the zenodo archive to download and extract.",
    type=click.STRING,
    default=None
)
@click.option(
    "-a",
    "--archive",
    help="The path to the extracted archive folder from zenodo.",
    type=click.Path(
        dir_okay=True, exists=True, path_type=pathlib.Path, file_okay=False
    ),
)
@click.option(
    "-p",
    "--partner-id",
    type=click.STRING,
    help="An identifier that should be applied to all of the extracted edges.",
)
@click.option(
    "-o",
    "--output",
    help="The name of the output folder the data should be wrote to.",
    type=click.Path(exists=False, path_type=pathlib.Path),
)
@click.option(
    "-e",
    "--experimental-data",
    help="The experimental csv file which should be used for public datasets.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
    default=None
)
def main(
    zenodo: str | None,
    archive: pathlib.Path,
    partner_id: str,
    output: pathlib.Path,
    experimental_data: pathlib.Path | None
):
    """
    Format the raw data extracted from the Zenodo archive into CSV files
    which can be used to process the results more easily.
    """
    # Make sure we are using pymabr3 for this script!
    version = pymbar.version.full_version.split(".")[0]
    if version != "3":
        raise RuntimeError(f"Pymbar3 should be used with this script, version: {pymbar.version.full_version} found")

    # create the output folder
    output.mkdir(parents=True, exist_ok=True)

    #  if we have an archive address download it
    if zenodo is not None:
        click.echo("Downloading archive from zenodo")
        # get the name of the archive folder
        archive_folder = download_zenodo_archive(address=zenodo, output_folder=output)
        click.echo("Extracting archive locally")
        shutil.unpack_archive(archive_folder, extract_dir=output)
        # get the name of the unpacked archive
        archive = archive_folder.with_suffix("")

    click.echo(f"Loading local archive {archive}")
    # assume unzipped archive has folders named after each dataset
    for dataset_name in archive.glob("*"):
        if dataset_name.suffix == ".zip":
            # unpack the results dir
            shutil.unpack_archive(dataset_name, extract_dir=dataset_name.with_suffix(""))
            # rename the dataset_name to match the unpacked archive
            dataset_name = dataset_name.with_suffix("")

        if dataset_name.is_dir():
            click.echo(f"Processing dataset {dataset_name}")
            # load the top level edge information from the all network properties file
            property_file = dataset_name.joinpath("all_network_properties.json")
            all_properties = json.loads(
                property_file.open("r").read(), cls=JSON_HANDLER.decoder
            )
            # find the experimental data csv file if not provided
            if experimental_data is not None:
                # assume it is a local FEP+ data file
                exp_df = pd.read_csv(experimental_data.as_posix())
                public_set = True
            else:
                exp_file =  list(dataset_name.glob("*.csv"))[0]
                # load the experimental data
                exp_df = load_exp_data(exp_file)
                public_set = False

            # extract the top level edge information into a df
            edge_data = extract_general_edge_data(
                all_properties=all_properties,
                dataset_name=dataset_name,
                industry_partner=partner_id,
                experimental_data=exp_df,
                archive=dataset_name,
                is_public=public_set
            )
            edge_data.to_csv(dataset_name.joinpath("pymbar3_edge_data.csv"))
            edge_data = pd.read_csv(dataset_name.joinpath("pymbar3_edge_data.csv"))
            # plot ddg vs exp
            click.echo("Plotting DDG vs Exp")
            plot_ddg_vs_experiment(
                edge_dataframe=edge_data,
                output=dataset_name,
                dataset_name=dataset_name.parent
            )
            click.echo("Plotting DDG between repeats")
            plot_ddg_repeats(
                edge_dataframe=edge_data,
                output=dataset_name
            )
            click.echo("Plotting DG vs Exp")
            dg_dataframe = calculate_and_plot_dg(
                edge_dataframe=edge_data,
                output=dataset_name,
                experimental_data=exp_df,
                dataset_name=dataset_name.parent,
                is_public=public_set
            )
            dg_dataframe.to_csv(dataset_name.joinpath("pymbar3_calculated_dg_data.csv"), index=False)
            # extract the cumulative estimates
            click.echo("Calculating cumulative data")
            cumulative_data = extract_cumulative_dg_estimates(
                all_properties=all_properties,
                dataset_name=dataset_name,
                industry_partner=partner_id,
                archive=dataset_name
            )
            cumulative_data.to_csv(dataset_name.joinpath("pymbar3_cumulative_data.csv"))


if __name__ == "__main__":
    main()
