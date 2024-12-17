# Dataset Analysis

The scripts provided in this folder should be used to process the industry benchmark datasets uploaded to Zenodo.

## Creating the analysis environments

You should use the recommended openfe-1.0.1 environment to run the analysis, this can be installed following the 
[guide](https://industrybenchmarks2024--157.org.readthedocs.build/en/157/installation.html#openfe-installation).

### Pymbar4
The Pymbar4 analysis script requires a seperate environment which can be created via:

```bash
mamba env craete -f pymbar4.yaml
```

## Running the analysis

The scripts should be executed in the following order.

``1_download_and_extract_data.py``: This script will download the results archive from Zenodo and extract edge scores 
and DG estimates for each repeat and phase and collect them into a CSV file (`pymbar3_edge_data.csv`) along with the experimental data. It will
also calculate the cumulative DG estimates for each edge and report these into a separate CSV file 
(`pymbar3_cumulative_data.csv`). 

The script will automatically find the experimental data for private sets, for public sets you can provide the data as:

```bash
python 1_download_and_extract_data.py -z https://zenodo.org/records/14229113   \ 
                                      -p Merck -o merck_hif2a                  \
                                      -e  hif2a_automap_symbmcorr_out.csv
```

`2_pymbar4.py`: This script will recalculate the DG, bootstrap error and overlap matrix for each edge in the extracted
archive using `pymbar4` and report the values into a new CSV file (`pymbar4_edge_data.csv`).




## Analysis Metadata

The following table should be used to check the meaning of each of the column names included in the CSV files. 

| Column Name                                     | Description                                                                                                              |
|-------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| `ligand_A/B`                                    | The identifier or name of the ligand                                                                                     |
| `lomap_score`                                   | The lomap score calculated for this transformation.                                                                      |
| `alchemical_charge_difference`                  | The difference in net formal charge between the ligands calculated as `ligand_A` - `ligand_B`.                           |
| `charge_score`                                  | A transformation difficulty score based on the net formal charge change.                                                 |
| `shape_score`                                   | A transformation difficulty score based on shape overlap of the ligands.                                                 |
| `volume_score`                                  | A transformation difficulty score based on the volume overlap of the ligands.                                            |
| `mapping_rmsd_score`                            | A transformation difficulty score based on the RMSD between the two ligands.                                             |
| `num_heavy_core`                                | The number of heavy atoms in the shared core topology between the two ligands.                                           |
| `num_heavy_dummy_A/B`                           | The number of unique heavy atoms in ligand A/B.                                                                          |
| `difference_num_rings_AB`                       | The difference in the number of rings between the ligands.                                                               |
| `difference_num_rot_bonds_AB`                   | The difference in the number of rotatable bonds between the ligands.                                                     |
| `morgan_tanimoto_similarity`                    | The tanimoto similarity of the Morgan fingerprints of the ligands.                                                       |
| `difference_solvent_accessible_surface_area`    | The difference in the solvent accessible surface area between the ligands.                                               |
| `atom_pair_dice_similarity`                     | The dice similarity of the atom pair fingerprints.                                                                       |
| `topological_torsion_dice_similarity`           | The dice similarity of the topological torsion fingerprints.                                                             |
| `dataset_name`                                  | The name of the dataset.                                                                                                 |
| `partner_id`                                    | The identifier of the industry partner that computed this dataset.                                                       |
| `exp DDG (kcal/mol)`                            | The experimental DDG in (kcal/mol)                                                                                       |
| `exp dDDG (kcal/mol)`                           | The experimental uncertainty in (kcal/mol)                                                                               |
| `solvent/complex_repeat_0/1/2_DG (kcal/mol)`    | The DG calculated for repeat (0/1/2) in phase (solvent/complex).                                                         |
| `solvent/complex_repeat_0/1/2_dDG (kcal/mol)`   | The error in the calculated DG for repeat (0/1/2) in phase (solvent/complex).                                            |
| `solvent/complex_repeat_0/1/2_com_drift_max`    | The maximum center of mass drift for the ligand in each repeat (0/1/2) and phase (solvent/complex).                      |
| `solvent/complex_repeat_0/1/2_ligand_rmsd_max`  | The maximum ligand RMSD in each repeat (0/1/2) and phase (solvent/complex).                                              |
| `solvent/complex_repeat_0/1/2_smallest_overlap` | The smallest off diagonal element of the mbar overlap matrix for each repeat (0/1/2) and phase (solvent/complex).        |
| `DG (kcal/mol)`                                 | The MLE calculated DG centered around 0 for the ligand.                                                                  |
| `uncertainty (kcal/mol)`                        | The uncertainty in the DG prediction for the ligand.                                                                     | 
| `phase`                                         | The phase of the transformation (for the cumulative DG values).                                                          |
| `repeat`                                        | The repeat of the transformation (for the cumulative DG values).                                                         |
| `Samples Xns DG`                                | The calculated DG for the given transformation leg using Xns of simulation data.                                         |
| `Samples Xns dDG`                               | The uncertainty in the DG prediction using Xns of simulation data.                                                       |
| `Samples Xns (subsample) DG`                    | The calculated DG for the transformation leg with subsampling calculated using Xns of simulation data.                   |
| `Samples Xns (subsample) dDG`                   | The uncertainty in the DG prediction with subsampling calculated using Xns of simulation data.                           |
| `Samples Xns (subsample_full) DG`               | The calculated DG for the transformation leg using Xns of simulation data with subsampling based on the full trajectory. |
| `Samples Xns (subsample_full) dDG`              | The uncertainty in the DG prediction using Xns of simulation data with subsampling based on the full trajectory.         |                                                                                            

