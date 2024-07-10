#!/usr/bin/env python

import argparse
import json
import os
import pathlib
import traceback
from pathlib import Path
from shutil import copyfile

import MDAnalysis as mda
import numpy as np
from openfe_analysis import FEReader
from openmmtools import multistate

# TODO test gather before and after


def compute_mbar_energies(analyzer):
    """
    Returns
    ------
    u_ln: energy matrix of shape (K,N) indexed by k,n
        K is the total number of states observables are desired.
        N is the total number of samples drawn from ALL states.
        The nth configuration is the energy evaluated in the kth thermodynamic state.
    N_l: 1-D iterable of shape K
        The number of samples drawn from each kth state.
    """
    analyzer.use_full_trajectory = True
    u_ln, N_l = analyzer._compute_mbar_decorrelated_energies()
    return u_ln, N_l


def get_replica_state_indices(analyzer):
    energy_data = list(analyzer._read_energies(truncate_max_n_iterations=True))
    replicas_state_indices = energy_data[-1]
    return replicas_state_indices


def subsample_traj(simulation, hybrid_system_pdb, lambda_windows, outfile):
    for i in range(0, lambda_windows):
        u = mda.Universe(hybrid_system_pdb, simulation, format=FEReader, state_id=i)
        frames = [round(i) for i in np.linspace(0, len(u.trajectory) - 1, 21)]
        out_traj = pathlib.Path(f"{outfile}_{i}.xtc")
        with mda.Writer(str(out_traj), n_atoms=len(u.atoms)) as w:
            for ts in u.trajectory[frames]:
                w.write(u.atoms)


def extract_data(simulation, checkpoint, hybrid_pdb, outfile, out_traj="out"):
    """
    Extract the MBAR-ready energy matrix and replica state indices.

    Parameters
    ----------
    simulation: pathlib.Path
        Path to the simulation `.nc` file
    checkpoint: pathlib.Path
        Path to the checkpoint `.chk` file
    hybrid_pdb: pathlib.Path
        Path to the `.pdb` file of the hybrid system
    outfile: pathlib.Path
        Path to the output `.npz` file
    out_traj: pathlib.Path
        Path to the output `.xtc` files. A separate file is created for every
        lambda window. The state number is appended to the filename. Default: 'out'
    """
    if not simulation.is_file() or not checkpoint.is_file():
        errmsg = "Either the simulation or checkpoint file could not be found"
        raise ValueError(errmsg)

    reporter = multistate.MultiStateReporter(
        storage=simulation.as_posix(),
        open_mode="r",
        checkpoint_storage=checkpoint.as_posix(),
    )

    analyzer = multistate.MultiStateSamplerAnalyzer(reporter)
    u_ln, N_l = compute_mbar_energies(analyzer)
    replicas_state_indices = get_replica_state_indices(analyzer)
    np.savez(outfile, u_ln=u_ln, N_l=N_l, replicas_state_indices=replicas_state_indices)

    lambda_windows = len(replicas_state_indices)
    subsample_traj(simulation, hybrid_pdb, lambda_windows, out_traj)


def make_backup(json_file: str) -> str:
    """
    Creates a backup of the given JSON file.
    The backup file is located next to the input file with BAK appended to the file name.

    Parameters
    ----------
    json_file : str
        Path to the JSON file to back up.

    Returns
    -------
    str
        Path to the backup file.
    """
    dir_name = os.path.dirname(json_file)
    base_name = os.path.basename(json_file)
    backup_file = os.path.join(dir_name, f"{base_name}.BAK")
    copyfile(json_file, backup_file)
    return backup_file


def restore_backup(json_file: str, backup_file: str) -> None:
    """
    Restores a JSON file from its backup.

    Parameters
    ----------
    json_file : str
        Path to the original JSON file.
    backup_file : str
        Path to the backup file.
    """
    copyfile(backup_file, json_file)


def delete_backup(backup_file: str) -> None:
    os.remove(backup_file)


def clean_results(json_files: list[str]) -> None:
    """
    Cleans up the results in the given JSON files.

    Note:
    To load in the structural analysis data

    >>> loaded = np.load("structural_analysis_data.npz")
    >>> loaded["time_ps"], loaded["protein_RMSD"], loaded["ligand_RMSD"], loaded["ligand_wander"], loaded["protein_2D_RMSD"]

    Parameters
    ----------
    json_files : List[str]
        List of paths to JSON files to clean up.
    """
    for json_file in json_files:
        if not os.path.exists(json_file):
            print(f"Error: {json_file} does not exist.")
            continue
        try:
            print(f"Working on file {json_file}")
            backup_file = make_backup(json_file)
            with open(json_file, "r") as f:
                results = json.load(f)

            # Check to see if we have already cleaned  up this result
            result_key = next(k for k in results["protocol_result"]["data"].keys())
            if (
                "structural_analysis"
                in results["protocol_result"]["data"][result_key][0]["outputs"]
            ):
                print("Cleaning up file")
            else:
                print("Skipping file, already cleaned")
                continue

            # Check to make sure we don't have more than one proto result
            # We might have ProtocolUnitResult-* and ProtocolUnitFailure-*
            # We only handel the case where we have one ProtocolUnitResult
            protocol_unit_result_count = len(
                [
                    k
                    for k in results["unit_results"].keys()
                    if k.startswith("ProtocolUnitResult")
                ]
            )

            # Check to make sure we don't just have failures
            # if all failures, tell user to re-run
            if protocol_unit_result_count != 1:
                print("More than one ProtocolUnitResult, skipping")
                continue
            elif protocol_unit_result_count == 0:
                print("All protocol units failed, skipping")
                continue

            # get the name of the key which is a gufe token
            # for the only ProtocolUnitResult-* in unit_results
            # this means we can grab the first that matches since there is only
            # one ProtocolUnitResult-*
            proto_key = next(
                k
                for k in results["unit_results"].keys()
                if k.startswith("ProtocolUnitResult")
            )
            results_dir = Path(
                results["unit_results"][proto_key]["outputs"]["nc"]["path"]
            ).parent
            structural_analysis_data = results["unit_results"][proto_key]["outputs"][
                "structural_analysis"
            ]

            # save structural analysis data
            np.savez_compressed(
                results_dir / "structural_analysis_data.npz",
                protein_RMSD=np.asarray(
                    structural_analysis_data["protein_RMSD"], dtype=np.float32
                ),
                ligand_RMSD=np.asarray(
                    structural_analysis_data["ligand_RMSD"], dtype=np.float32
                ),
                ligand_wander=np.asarray(
                    structural_analysis_data["ligand_wander"], dtype=np.float32
                ),
                protein_2D_RMSD=np.asarray(
                    structural_analysis_data["protein_2D_RMSD"], dtype=np.float32
                ),
                time_ps=np.asarray(
                    structural_analysis_data["time(ps)"], dtype=np.int32
                ),
            )

            # remove structural_analysis data stuffed into unit results
            del results["unit_results"][proto_key]["outputs"]["structural_analysis"]

            # Now we subsamble the traj and save reporter data
            simulation = results_dir / "simulation.nc"
            checkpoint = results_dir / "checkpoint.chk"
            hybrid_pdb = results_dir / "hybrid_system.pdb"
            # TODO better name?
            outfile = results_dir / "multi_state_sampler_data.npz"
            out_traj = results_dir / "out"
            extract_data(simulation, checkpoint, hybrid_pdb, outfile, out_traj="out")

            # Now we delete files we don't need anymore
            os.remove(simulation)
            os.remove(checkpoint)

            # remove structural_analysis data stuffed into protocol_result
            # and remove ligand + pdb
            # this data is duped in structural_analysis_data.npz
            del results["protocol_result"]["data"][result_key][0]["outputs"][
                "structural_analysis"
            ]
            # do we want to keep these or not? lig and pdb could be useful, but we have the inputs
            # I think Irfan said keep
            del results["protocol_result"]["data"][result_key][0]["inputs"]["stateA"]
            del results["protocol_result"]["data"][result_key][0]["inputs"]["stateB"]
            del results["protocol_result"]["data"][result_key][0]["inputs"][
                "ligandmapping"
            ]

            # TODO save as gzip -- maybe, gather will fail then?
            with open(json_file, "w") as f:
                json.dump(results, f)

        except Exception as e:
            print("oh no, we hit an error, restoring backup")
            restore_backup(json_file, backup_file)
            print(traceback.format_exc())
            raise e

        else:
            print("removing backup")
            delete_backup(backup_file)


def main():
    parser = argparse.ArgumentParser(
        description="Cleanup results directory and results json"
    )
    parser.add_argument(
        "json_files", metavar="JSON_FILE", nargs="+", help="JSON file(s) to reduce size"
    )

    args = parser.parse_args()

    clean_results(args.json_files)


if __name__ == "__main__":
    main()
