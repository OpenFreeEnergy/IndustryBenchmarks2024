#!/usr/bin/env python

import json
import numpy as np
import os
from shutil import copyfile
import argparse
import traceback

# the plan
# read in result json to get path to results dir
# remove files we don't need in results dir
# trim down the traj file in results dir
# make a backup of results json
# trim down results json and overwerite the old one
# this means splitting out the data into numpy arrays
# delete the backup (use try, except, finally)

# TODO get this from argparse
json_file = "results_tyk2_transformations_openfe_1.0.0rc0_uc_0/lig_jmc_23_complex_lig_jmc_28_complex_complex.json"

def make_backup(json_file):
    dir_name = os.path.dirname(json_file)
    base_name = os.path.basename(json_file)
    backup_file = os.path.join(dir_name, f"{base_name}.BAK")
    copyfile(json_file, backup_file)
    return backup_file

def restore_backup(json_file, backup_file):
    copyfile(backup_file, json_file)

def delete_backup(backup_file):
    os.remove(backup_file)

try:
    backup_file = make_backup(json_file)
    with open(json_file, 'r') as f:
        results = json.load(f)
    # remove pdb + ligand stuffed into the result
    del results["protocol_result"]["data"]
    # Check to make sure we don't have more than one proto result
    assert len(results["unit_results"].keys()) == 1
    # get the name of the key which is a gufe token
    proto_key = next(iter(results["unit_results"]))
    keys = results["unit_results"][proto_key]["outputs"]["structural_analysis"]
    structural_analysis_data = results["unit_results"][proto_key]["outputs"]["structural_analysis"]
    # get this data back like
    # loaded = np.load("structural_analysis_data.npz")
    # loaded["time_ps"], loaded["protein_RMSD"], loaded["ligand_RMSD"], loaded["ligand_wander"], loaded["protein_2D_RMSD"]
    # TODO save this in the results directory, which we can figure out from the results json
    np.savez_compressed("structural_analysis_data.npz",
                    protein_RMSD=structural_analysis_data["protein_RMSD"],
                    ligand_RMSD=structural_analysis_data["ligand_RMSD"],
                    ligand_wander=structural_analysis_data["ligand_wander"],
                    protein_2D_RMSD=structural_analysis_data["protein_2D_RMSD"],
                    time_ps=structural_analysis_data["time(ps)"]
                   )
    del results["unit_results"][proto_key]["outputs"]["structural_analysis"]
    with open(json_file, 'w') as f:
        json.dump(results, f)
except Exception as e:
    print("oh no, we hit an error, restoring backup")
    restore_backup(json_file, backup_file)
    print(traceback.format_exc())
    raise e
else:
    print("removing backup")
    delete_backup(backup_file)
