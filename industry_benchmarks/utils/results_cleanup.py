#!/usr/bin/env python

import argparse
import json
import os
import traceback
from pathlib import Path
from shutil import copyfile

# TODO merge into a single file if not used elsewhere
from traj_cleanup import extract_data

import numpy as np

# TODO
# remove files we don't need in results dir
# trim down the traj file in results dir
# use traj_cleanup.py to sub sample, then delete nc and chk
# only delete out this bit from data
# d["protocol_result"]["data"]["135312977195331644756224466752878866009"][0]["outputs"]["structural_analysis"]

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
            if "data" in results["protocol_result"]:
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
                print("All protocol units failed")
                continue
            # get the name of the key which is a gufe token
            # TODO make sure we don't grab a result faiulre
            proto_key = next(iter(results["unit_results"]))
            results_dir = Path(
                results["unit_results"][proto_key]["outputs"]["nc"]["path"]
            ).parent
            structural_analysis_data = results["unit_results"][proto_key]["outputs"][
                "structural_analysis"
            ]
            # save structural analysis data
            # TODO save as 32 bit floats
            np.savez_compressed(
                results_dir / "structural_analysis_data.npz",
                protein_RMSD=structural_analysis_data["protein_RMSD"],
                ligand_RMSD=structural_analysis_data["ligand_RMSD"],
                ligand_wander=structural_analysis_data["ligand_wander"],
                protein_2D_RMSD=structural_analysis_data["protein_2D_RMSD"],
                time_ps=structural_analysis_data["time(ps)"],
            )
            # remove structural_analysis data stuffed into results
            del results["unit_results"][proto_key]["outputs"]["structural_analysis"]
            for key in results["protocol_result"]["data"]:
                del results["protocol_result"]["data"][key][0]["outputs"]["structural_analysis"]
            # remove pdb + ligand stuffed into the result
            # TODO remove
            # del results["protocol_result"]["data"]
            # TODO save as gzip
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
