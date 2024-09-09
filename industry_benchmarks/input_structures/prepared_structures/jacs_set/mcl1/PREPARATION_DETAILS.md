# GSK Industry Benchmark
Prepared by Alexander Williams
## MCL1 Set
## Software Used
1. Maestro Protein Prep v.2024-02
## Protein Preparation: Steps Taken
1. Protein from original structures was imported into Maestro and using the protein preparation panel the interactive mode was activated following settings were applied on the preprocess panel.
### Maestro Protein Preparation Workflow Selected Settings
   1. Cap Termini
   2. Convert selenomethionines to methionines
   3. Include Peptides when capping termini.
2. Final NMA peptide prepared by maestro was converted to NME by using the builder panel and selecting other edits ... change atom properties.
3. Protein was saved out to protein.pdb

## Ligand Preparation: Steps Taken
1. Utilizing the included edges csv file within the original_structures folder, each ligand found within the included *ligands.sdf file was confirmed to have a node within said edges file. These ligands were then loaded into maestro to confirm their structure and saved to ligands.sdf.
