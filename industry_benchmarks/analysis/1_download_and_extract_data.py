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
from concurrent.futures import ProcessPoolExecutor, as_completed
from openmmtools.multistate import MultiStateSamplerAnalyzer, MultiStateReporter, utils
import logging
import tempfile


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
    # we need to replace the -1 values which indicate the error was not supplied with 0.0
    exp_data["Affinity Error (nM)"] = exp_data["Affinity Error (nM)"].mask(exp_data["Affinity Error (nM)"] == -1, 0.0)
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
    # make it clear which result can not be found!
    try:
        ligand_a_dg = experimental_data[
            experimental_data["Ligand name"] == ligand_a
            ].iloc[0]["Exp. dG (kcal/mol)"]
        if "Exp. dG error (kcal/mol)" in experimental_data.columns:
            ligand_a_ddg = experimental_data[
                experimental_data["Ligand name"] == ligand_a
                ].iloc[0]["Exp. dG error (kcal/mol)"]
        else:
            ligand_a_ddg = 0.0
    except IndexError as e:
        print(ligand_a, " not found!")
        raise e
    try:
        ligand_b_dg = experimental_data[
            experimental_data["Ligand name"] == ligand_b
            ].iloc[0]["Exp. dG (kcal/mol)"]
        if "Exp. dG error (kcal/mol)" in experimental_data.columns:
            ligand_b_ddg = experimental_data[
                experimental_data["Ligand name"] == ligand_b
            ].iloc[0]["Exp. dG error (kcal/mol)"]
        else:
            ligand_b_ddg = 0.0
    except IndexError as e:
        print(ligand_b, " not found!")
        raise e
    ddg = ligand_b_dg - ligand_a_dg
    ddg_error = (ligand_a_ddg ** 2 + ligand_b_ddg ** 2) ** 0.5
    return ddg, ddg_error


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

def _extract_edge_data(
        edge_data: dict[str, float],
        all_props_data: dict,
        archive: pathlib.Path,
) -> dict:
    for phase in ["solvent", "complex"]:
        for repeat in range(3):
            repeat_name = f"{phase}_{edge_data['ligand_A']}_{edge_data['ligand_B']}_repeat_{repeat}"
            phase_repeat_name = f"{phase}_repeat_{repeat}"
            try:
                dg, error = all_props_data["DDG_estimates"][repeat_name]
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
            except KeyError:
                # sometimes we have edges that are not simulated see merck cmet
                continue
    return edge_data


def extract_general_edge_data(
    all_properties: dict[str, dict],
    dataset_name: str,
    industry_partner: str,
    experimental_data: pd.DataFrame,
    archive: pathlib.Path,
    is_public: bool,
    workers: int
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
    workers: int
        The number of workers to use when extracting the edge data.

    Returns
    -------
        A DataFrame of the edge data, such as edge scores and DG values per phase and repeat.

    Notes
    -----
    Running MBAR a lot of times is slow so use more workers where possible.
    """
    all_data = []
    job_list = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for edge_data in all_properties["Edges"].values():
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
            job_list.append(pool.submit(_extract_edge_data, edge_data, all_properties, archive))

        # unpack the work
        for job in tqdm.tqdm(as_completed(job_list), desc="Extracting edge data", total=len(job_list)):
            result = job.result()
            if "solvent_repeat_0_DG (kcal/mol)" in result:
                # make sure there was data for this edge
                all_data.append(result)

    return pd.DataFrame(all_data)


def download_zenodo_archive(address: str, output_folder: pathlib.Path) -> list[str]:
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
    The names of the files downloaded from the archive
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
            files = []
            for f in record_files:
                filename = f.get("filename") or f["key"]
                file_address = f"https://zenodo.org/records/{record_id}/files/{filename}"
                wget.download(file_address, output_folder.as_posix())
                files.append(filename)
            return files
        else:
            record_data.raise_for_status()

    else:
        raise RuntimeError(f"The zenodo address ({address}) did not contain the record id")



def get_edge_matrix(
    edge_name: str,
    archive: pathlib.Path,
    matrix_type: Literal[
        "com", "ligand_rmsd", "protein_rmsd", "reduced_potential", "samples", "state_index"
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
        "com": ("ligand_wander.txt", np.float64),
        "ligand_rmsd": ("ligand_RMSD.txt", np.float64),
        "protein_rmsd": ("protein_RMSD.txt", np.float64),
        "reduced_potential": ("u_ln.txt", np.float64),
        "samples": ("N_l.txt", np.int16),
        "state_index": ("replicas_state_indices.txt", np.int16)
    }
    file_name = archive.joinpath(edge_name, types_to_files[matrix_type][0])
    return np.loadtxt(file_name).astype(types_to_files[matrix_type][1])


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
        samples: np.array,
        state_indices: np.array,
        charge_change: bool,

) -> dict[str, float]:
    """
    Compute the cumulative DG estimates for each nanosecond of simulation time with and without subsampling.

    Also calculate the DG for each nanosecond using the full simulation autocorrelation time.

    Parameters
    ----------
    reduced_potential: np.array
        The reduced potential calculated for each sample in the form (n_states, n_samples).
    samples: np.array
        The number of samples per state.
    state_indices: np.array
        The index of the state each sample is calculated at.
    charge_change: bool
        If the edge involves a charge change which changes the amount of simulation time and lambda windows.

    Returns
    -------
    A dict of the nanoseconds of simulation time and the estimated DG and dDG in kcal/mol.
    """
    # hide the openmmtools warnings
    logger = logging.getLogger()
    logger.setLevel(logging.ERROR)

    # fake the reporter and analyzer
    temp_file = tempfile.NamedTemporaryFile().name
    reporter = MultiStateReporter(open_mode="w", storage=temp_file)
    analyzer = MultiStateSamplerAnalyzer(reporter)

    cumulative_dgs = {}

    # assuming we ran at the default openfe temperature in the protocol
    kt = unit.molar_gas_constant * 298.15 * unit.kelvin

    # workout how many nanoseconds were run
    total_nano = 5 if not charge_change else 20
    # get the number of replicas run
    n_replicas = len(samples)
    # get the number of samples from each state
    n_states = int(samples[0])

    # make sure we have the right number of replicas
    if total_nano == 5:
        assert n_replicas == 11
    else:
        assert n_replicas == 22

    # format the reduced potential to work with openmmtools
    energy_matrix = np.zeros((n_replicas, n_replicas, n_states))
    for i in range(n_replicas):
        energy_matrix[i] = reduced_potential[:, i * n_states: int(i + 1) * n_states]

    # calculate the statistical inefficiency for the full simulation
    analyzer.clear()
    full_sim_number_equilibrated, full_sim_g_t, _ = analyzer._get_equilibration_data(energy_matrix,
                                                                          replica_state_indices=state_indices)

    # calculate the DG for each nanosecond
    for i in range(1, total_nano + 1):
        analyzer.clear()
        # get the % of data to slice, this matches the simulation_real_time_analysis.yaml exactly
        n_samples = int(np.ceil(n_states * (i / total_nano)))
        # slice out the data we would have at this simulation time
        sample_matrix = energy_matrix[:,:,:n_samples]
        sample_state_matrix = state_indices[:, :n_samples]
        # calculate dg without subsamples
        sample_u_ln_matrix = analyzer.reformat_energies_for_mbar(sample_matrix)
        samples_per_state = np.array([n_samples for _ in range(n_replicas)], dtype=int)

        mbar = MBAR(sample_u_ln_matrix, samples_per_state)
        try:
            estimate = mbar.getFreeEnergyDifferences(return_dict=True)
            DG, dDG = estimate["Delta_f"][0, -1] * kt, estimate["dDelta_f"][0, -1] * kt
            cumulative_dgs[f"Samples {i}ns DG"] = DG.to(unit.kilocalorie_per_mole).m
            cumulative_dgs[f"Samples {i}ns dDG"] = dDG.to(unit.kilocalorie_per_mole).m
        except pymbar.utils.ParameterError:
            # this is raised if we have an error solving
            # in this case just return na?
            cumulative_dgs[f"Samples {int(i * 10)}% DG"] = pd.NA
            cumulative_dgs[f"Samples {int(i * 10)}% dDG"] = pd.NA

        # again with subsample for partial data
        analyzer.clear()
        # use openmmtools to subsample the data
        number_equilibrated, g_t, _ = analyzer._get_equilibration_data(sample_matrix, replica_state_indices=sample_state_matrix)
        # remove unequilibrated data
        sample_matrix = utils.remove_unequilibrated_data(sample_matrix, number_equilibrated, -1)
        sample_state_matrix = utils.remove_unequilibrated_data(sample_state_matrix, number_equilibrated, -1)
        # subsample using g_t
        sample_matrix = utils.subsample_data_along_axis(sample_matrix, g_t, -1)
        sample_state_matrix = utils.subsample_data_along_axis(sample_state_matrix, g_t, -1)

        # Determine how many samples and which states they were drawn from.
        unique_sampled_states, counts = np.unique(sample_state_matrix, return_counts=True)
        # Assign those counts to the correct range of states.
        samples_per_state = np.zeros([n_replicas], dtype=int)
        samples_per_state[:][unique_sampled_states] = counts
        sample_u_ln_matrix = analyzer.reformat_energies_for_mbar(sample_matrix)
        mbar = MBAR(sample_u_ln_matrix, samples_per_state)

        # run mbar again
        try:
            estimate = mbar.getFreeEnergyDifferences(return_dict=True)
            DG, dDG = estimate["Delta_f"][0, -1] * kt, estimate["dDelta_f"][0, -1] * kt
            cumulative_dgs[f"Samples {i}ns (subsample) DG"] = DG.to(unit.kilocalorie_per_mole).m
            cumulative_dgs[f"Samples {i}ns (subsample) dDG"] = dDG.to(unit.kilocalorie_per_mole).m
        except pymbar.utils.ParameterError:
            # this is raised if we can not solve
            cumulative_dgs[f"Samples {i}ns (subsample) DG"] = pd.NA
            cumulative_dgs[f"Samples {i}ns (subsample) dDG"] = pd.NA

        # again with subsample with full g_t
        analyzer.clear()
        # get the % of data to slice, this matches the simulation_real_time_analysis.yaml exactly
        n_samples = int(np.ceil(n_states * (i / total_nano)))
        # slice out the data we would have at this simulation time
        sample_matrix = energy_matrix[:, :, :n_samples]
        sample_state_matrix = state_indices[:, :n_samples]
        # remove unequilibrated data using the full simulation g_t
        sample_matrix = utils.remove_unequilibrated_data(sample_matrix, full_sim_number_equilibrated, -1)
        sample_state_matrix = utils.remove_unequilibrated_data(sample_state_matrix, full_sim_number_equilibrated, -1)
        # subsample using g_t
        sample_matrix = utils.subsample_data_along_axis(sample_matrix, full_sim_g_t, -1)
        sample_state_matrix = utils.subsample_data_along_axis(sample_state_matrix, full_sim_g_t, -1)

        # Determine how many samples and which states they were drawn from.
        unique_sampled_states, counts = np.unique(sample_state_matrix, return_counts=True)
        # Assign those counts to the correct range of states.
        samples_per_state = np.zeros([n_replicas], dtype=int)
        samples_per_state[:][unique_sampled_states] = counts
        sample_u_ln_matrix = analyzer.reformat_energies_for_mbar(sample_matrix)

        # run mbar again
        try:
            mbar = MBAR(sample_u_ln_matrix, samples_per_state)
            estimate = mbar.getFreeEnergyDifferences(return_dict=True)
            DG, dDG = estimate["Delta_f"][0, -1] * kt, estimate["dDelta_f"][0, -1] * kt
            cumulative_dgs[f"Samples {i}ns (subsample_full) DG"] = DG.to(unit.kilocalorie_per_mole).m
            cumulative_dgs[f"Samples {i}ns (subsample_full) dDG"] = dDG.to(unit.kilocalorie_per_mole).m
        except (pymbar.utils.ParameterError, IndexError):
            # this is raised if we can not solve
            # or if we have removed all samples from the matrix
            cumulative_dgs[f"Samples {i}ns (subsample_full) DG"] = pd.NA
            cumulative_dgs[f"Samples {i}ns (subsample_full) dDG"] = pd.NA

    return cumulative_dgs

def extract_cumulative_dg_estimates(
    all_properties: dict[str, dict],
    industry_partner: str,
    archive: pathlib.Path,
    workers: int
) -> pd.DataFrame:
    """
    For each edge in the dataset estimate the DG cumulative values after each nanosecond of simulation time.

    Parameters
    ----------
     all_properties: dict
        A dictionary of all network properties created by the data_gathering.py script.
    industry_partner: str
        The unique name which should be used to identify the partner this dataset came from.
    archive: pathlib.Path
        The location of the archive to extract edge matrix data from.
    workers: int
        The number of workers to use when extracting the edge data.

    Returns
    -------
        A dataframe of DG estimates using portions of the total data, where each row is a unique phase and repeat.

    Notes
    -----
    Running MBAR a lot of times is very expensive use as many workers as possible.
    We assume 5 nanoseconds and 11 windows for no charge change and 22 windows and 20 nanoseconds for a charge change.
    """
    all_data = []
    job_list = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for edge_data in all_properties["Edges"].values():
            # add in the per repeat value for each phase
            job_list.append(pool.submit(_get_cumulative_edge_estimate, edge_data, archive, industry_partner))

        for result in tqdm.tqdm(as_completed(job_list), desc="Calculating cumulative DGs", total=len(job_list)):
            result = result.result()
            if result:
                all_data.extend(result)

    return pd.DataFrame(all_data)

def _get_cumulative_edge_estimate(
        edge_data: dict,
        archive: pathlib.Path,
        partner_id: str
) -> list[dict[str, float]]:

    cumulative_data = []
    for phase in ["solvent", "complex"]:
        for repeat in range(3):
            repeat_name = f"{phase}_{edge_data['ligand_A']}_{edge_data['ligand_B']}_repeat_{repeat}"
            # collect some general data
            try:
                dg_data = {
                    "ligand_A": edge_data["ligand_A"],
                    "ligand_B": edge_data["ligand_B"],
                    "phase": phase,
                    "repeat": repeat,
                    "partner_id": partner_id
                }
                charge_change = False if edge_data["charge_score"] == 1.0 else True
                potential = get_edge_matrix(
                    edge_name=repeat_name,
                    archive=archive,
                    matrix_type="reduced_potential",
                )
                samples = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="samples"
                )
                state_indices = get_edge_matrix(
                    edge_name=repeat_name, archive=archive, matrix_type="state_index"
                )
                cumulative_dg = get_cumulative_dg(
                    reduced_potential=potential.T,
                    samples=samples,
                    charge_change=charge_change,
                    state_indices=state_indices.T
                )
                # merge the cumulative data
                dg_data.update(cumulative_dg)
                # save the data for this transformation
                cumulative_data.append(dg_data)
            except FileNotFoundError:
                # some transforms have no results
                continue
    return cumulative_data


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
    try:
        plot_DGs(
            graph=fe_map.to_legacy_graph(),
            title=dataset_name,
            filename=output.joinpath("DG_vs_exp.svg").as_posix(),
            figsize=4
        )
    except ValueError:
        # if to few results to get a linear correlation plot without it
        temp_graph = fe_map.to_legacy_graph()
        x_data = np.asarray([node[1]["exp_DG"] for node in temp_graph.nodes(data=True)])
        y_data = np.asarray([node[1]["calc_DG"] for node in temp_graph.nodes(data=True)])
        xerr = np.asarray([node[1]["exp_dDG"] for node in temp_graph.nodes(data=True)])
        yerr = np.asarray([node[1]["calc_dDG"] for node in temp_graph.nodes(data=True)])
        # centralise
        x_data = x_data - np.mean(x_data)
        y_data = y_data - np.mean(y_data)
        _master_plot(
            x_data,
            y_data,
            xerr=xerr,
            yerr=yerr,
            origins=False,
            statistics=["RMSE", "MUE", "rho"],
            quantity=rf"$\Delta$ G",
            title=dataset_name,
            filename=output.joinpath("DG_vs_exp.svg").as_posix(),
            figsize=4,
            bootstrap_x_uncertainty=False,
            bootstrap_y_uncertainty=False,
            statistic_type="mle",
        )
    # if we have a public result plot DG vs FEP+
    absolute_df = fe_map.get_absolute_dataframe()
    abs_calc = absolute_df[absolute_df["computational"] == True].copy(deep=True)
    exp_data = absolute_df[absolute_df["computational"] == False]

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

    # make a new dataframe with the predicted and experimental DGs shifted by the mean of the experimental data
    # shift the predicted results
    mean_exp = exp_data["DG (kcal/mol)"].mean()
    abs_calc["DG (kcal/mol)"] += mean_exp
    # loop over and add the exp data
    exp_values, exp_error = [], []
    for _, row in abs_calc.iterrows():
        exp_row = exp_data[exp_data["label"] == row["label"]].iloc[0]
        exp_values.append(exp_row["DG (kcal/mol)"])
        exp_error.append(exp_row["uncertainty (kcal/mol)"])
    abs_calc["Exp DG (kcal/mol)"] = exp_values
    abs_calc["Exp dDG (kcal/mol)"] = exp_error
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
@click.option(
    "-w",
    "--workers",
    default=2,
    type=click.INT,
    help="The number of workers to use to extract the edge data."
)
@click.option(
    "--download-only",
    default=False,
    is_flag=True,
    help="If the Zenodo archive should be downloaded and not processed, needed when multiple results are in a single archive. "
)
def main(
    zenodo: str | None,
    archive: pathlib.Path,
    partner_id: str,
    output: pathlib.Path,
    experimental_data: pathlib.Path | None,
    workers: int,
    download_only: bool
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
        archive_files = download_zenodo_archive(address=zenodo, output_folder=output)
        # unpack each one
        for file in archive_files:
            shutil.unpack_archive(output.joinpath(file), extract_dir=output)

        # Script can only do one archive at a time due to exp data, so quit if more than one
        if len(archive_files) > 1:
            exit()
        else:
            archive = output.joinpath(archive_files[0]).with_suffix("")

    if download_only:
        exit()

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
                exp_df = pd.read_csv(experimental_data.as_posix(), dtype={"Ligand name": str})
                public_set = True
            else:
                exp_file =  list(dataset_name.parent.glob("*.csv"))[0]
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
                is_public=public_set,
                workers=workers
            )
            edge_data.to_csv(dataset_name.joinpath("pymbar3_edge_data.csv"))
            # grab the non-failed edges if we have failures
            if "failed" in edge_data.columns:
                complete_df = edge_data[edge_data["failed"] != True].copy(deep=True)
                complete_df.reset_index(inplace=True)
            else:
                complete_df = edge_data.copy(deep=True)
            # plot ddg vs exp
            click.echo("Plotting DDG vs Exp")
            plot_ddg_vs_experiment(
                edge_dataframe=complete_df,
                output=dataset_name,
                dataset_name=dataset_name.parent
            )
            click.echo("Plotting DDG between repeats")
            plot_ddg_repeats(
                edge_dataframe=complete_df,
                output=dataset_name
            )
            click.echo("Plotting DG vs Exp")
            dg_dataframe = calculate_and_plot_dg(
                edge_dataframe=complete_df,
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
                industry_partner=partner_id,
                archive=dataset_name,
                workers=workers
            )
            cumulative_data.to_csv(dataset_name.joinpath("pymbar3_cumulative_data.csv"))


if __name__ == "__main__":
    main()

