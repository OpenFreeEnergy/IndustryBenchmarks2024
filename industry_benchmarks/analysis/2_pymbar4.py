import click
import pathlib
import pandas as pd
import pymbar
import tqdm
from pymbar import MBAR
import json
import numpy as np
from gufe.tokenization import JSON_HANDLER
from typing import Literal
from openff.units import unit
import logging
import tempfile
from openmmtools.multistate import MultiStateSamplerAnalyzer, MultiStateReporter, utils
import yaml



def compute_mabr_values(
        reduced_potential: np.array,
        samples: np.array,
        n_equilibration_iterations: int,
        statistical_inefficiency: float,
        state_indices: np.array
) -> dict[str, np.array]:
    """
    Compute the MBAR overlap matrix using pymbar4 as well as DG and the bootstrap error.

    Parameters
    ----------
    reduced_potential: np.array
        The reduced potential calculated for each sample in the form (n_states, n_samples)
    samples
        The number of samples per state.
    n_equilibration_iterations: int
        The number of equilibrium iteration to discard
    statistical_inefficiency: float
        The statistical_inefficiency of the production run used to subsample the trajectory
    state_indices: np.array
        The index of the state each sample is calculated at
    Returns
    -------
        The MBAR overlap matrix calculated with pymbar.
    """
    # hide the openmmtools warnings
    logger = logging.getLogger()
    logger.setLevel(logging.ERROR)

    # fake the reporter and analyzer
    temp_file = tempfile.NamedTemporaryFile().name
    reporter = MultiStateReporter(open_mode="w", storage=temp_file)
    analyzer = MultiStateSamplerAnalyzer(reporter)
    analyzer.clear()

    # format the reduced potential to work with openmmtools
    n_replicas = len(samples)
    n_states = int(samples[0])
    energy_matrix = np.zeros((n_replicas, n_replicas, n_states))
    for i in range(n_replicas):
        energy_matrix[i] = reduced_potential[:, i * n_states: int(i + 1) * n_states]

    # remove unequilibrated data
    sample_matrix = utils.remove_unequilibrated_data(energy_matrix, n_equilibration_iterations, -1)
    sample_state_matrix = utils.remove_unequilibrated_data(state_indices, n_equilibration_iterations, -1)
    # subsample using g_t
    sample_matrix = utils.subsample_data_along_axis(sample_matrix, statistical_inefficiency, -1)
    sample_state_matrix = utils.subsample_data_along_axis(sample_state_matrix, statistical_inefficiency, -1)

    # Determine how many samples and which states they were drawn from.
    unique_sampled_states, counts = np.unique(sample_state_matrix, return_counts=True)
    # Assign those counts to the correct range of states.
    samples_per_state = np.zeros([n_replicas], dtype=int)
    samples_per_state[:][unique_sampled_states] = counts
    sample_u_ln_matrix = analyzer.reformat_energies_for_mbar(sample_matrix)

    mbar = MBAR(sample_u_ln_matrix, samples_per_state, solver_protocol="robust", n_bootstraps=1000, bootstrap_solver_protocol="robust")
    dg_data = mbar.compute_free_energy_differences(uncertainty_method="bootstrap", compute_uncertainty=True)
    overlap_matrix = mbar.compute_overlap()
    dg_data["overlap"] = get_smallest_off_diagonal(overlap_matrix["matrix"])
    return dg_data

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

def get_edge_matrix(
    edge_name: str,
    archive: pathlib.Path,
    matrix_type: Literal[
        "reduced_potential", "samples", "state_index"
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
        "reduced_potential": ("u_ln.txt", np.float64),
        "samples": ("N_l.txt", np.int16),
        "state_index": ("replicas_state_indices.txt", np.int16)
    }
    file_name = archive.joinpath(edge_name, types_to_files[matrix_type][0])
    return np.loadtxt(file_name).astype(types_to_files[matrix_type][1])

def get_subsample_data(edge_name: str, archive: pathlib.Path) -> tuple[int, float]:
    """
    Get the calculated number of equilibrium iterations to remove and the
    statistical_inefficiency for this transformation.

    Parameters
    ----------
    edge_name: str
        The name of the edge of the format phase_ligandA_ligand_B_repeat_x
    archive: pathlib.Path
        The path to the extracted zenodo archive.

    Returns
    -------
    A tuple of the number of iterations to remove and the statistical_inefficiency.
    """
    data = yaml.safe_load(archive.joinpath(edge_name, "simulation_real_time_analysis.yaml").read_text())
    return data[-1]["mbar_analysis"]["n_equilibrium_iterations"], data[-1]["mbar_analysis"]["statistical_inefficiency"]


def get_edge_data(edge_data: dict, dataset_path: pathlib.Path) -> dict[str, float]:
    """
    For the given edge compute the DG, bootstrap error and the smallest overlap matrix element for each repeat and phase

    Parameters
    ----------
    edge_data: dict
        The edge data from the all properties dict
    dataset_path: pathlib.Path
        The path to this dataset that we should extract the reduced potential from.

    Returns
    -------
    A dict DG, bootstrap error and smallest overlap matrix element per phase and repeat.
    """
    kt = unit.molar_gas_constant * 298.15 * unit.kelvin
    return_data = {
        "ligand_A": edge_data["ligand_A"],
        "ligand_B": edge_data["ligand_B"]
    }
    for phase in ["solvent", "complex"]:
        for repeat in range(3):
            try:
                repeat_name = f"{phase}_{edge_data['ligand_A']}_{edge_data['ligand_B']}_repeat_{repeat}"
                phase_repeat_name = f"{phase}_repeat_{repeat}"
                potential = get_edge_matrix(
                    edge_name=repeat_name,
                    archive=dataset_path,
                    matrix_type="reduced_potential",
                )
                samples = get_edge_matrix(
                    edge_name=repeat_name, archive=dataset_path, matrix_type="samples"
                )
                state_indices = get_edge_matrix(
                    edge_name=repeat_name, archive=dataset_path, matrix_type="state_index"
                )
                n_equil, statistical_inefficiency = get_subsample_data(
                    edge_name=repeat_name, archive=dataset_path
                )
                dg_data = compute_mabr_values(
                    reduced_potential=potential.T,
                    samples=samples,
                    state_indices=state_indices.T,
                    n_equilibration_iterations=n_equil,
                    statistical_inefficiency=statistical_inefficiency
                )
                # extract the DG, bootstrap error and smallest overlap
                DG, dDG = dg_data["Delta_f"][0, -1] * kt, dg_data["dDelta_f"][0, -1] * kt
                return_data[f"{phase_repeat_name}_DG (kcal/mol)"] = DG.to(unit.kilocalorie_per_mole).m
                return_data[f"{phase_repeat_name}_dDG (kcal/mol)"] = dDG.to(unit.kilocalorie_per_mole).m
                return_data[f"{phase_repeat_name}_smallest_overlap"] = dg_data["overlap"]
            except FileNotFoundError:
                continue

    return return_data

@click.command()
@click.option(
    "-a",
    "--archive",
    help="The path to the extracted archive folder from zenodo.",
    type=click.Path(
        dir_okay=True, exists=True, path_type=pathlib.Path, file_okay=False)
)
def main(archive: pathlib.Path):
    """
    Run Pymbar4 and calculate DG the bootstrapped error and the overlap matrix for the set of edges and write out the
    results to a CSV in the same folder.
    """
    version = pymbar.__version__.split(".")[0]
    if version != "4":
        raise RuntimeError(f"Pymbar4 should be used with this script, version: {pymbar.__version__} found")

    click.echo(f"Loading local archive {archive}")
    # assume unzipped archive has folders named after each dataset
    for dataset_name in archive.glob("*"):
        if dataset_name.is_dir():
            click.echo(f"Processing dataset {dataset_name}")
            # load the top level edge information from the all network properties file
            property_file = dataset_name.joinpath("all_network_properties.json")
            all_properties = json.loads(
                property_file.open("r").read(), cls=JSON_HANDLER.decoder
            )

            all_edge_data = []
            for edge_data in tqdm.tqdm(all_properties["Edges"].values(), desc="Calculating pymbar4", total=len(all_properties["Edges"])):
                result = get_edge_data(edge_data, dataset_name)
                if len(result) == 20:
                    all_edge_data.append(result)
                # if not some results are missing so skip

            df = pd.DataFrame(all_edge_data)

            df.to_csv(dataset_name.joinpath("pymbar4_edge_data.csv"))



if __name__ == "__main__":
    main()