# GSK Industry Benchmark
Prepared by Alexander Williams
## HNE MCS Set
## Software Used
1. Maestro Protein Prep v.2024-02
2. Cresset Flare v8
## Protein Preparation: Steps Taken
1. Protein from original structures was imported into Maestro and using the protein preparation panel the interactive mode was activated following settings were applied on the preprocess panel.
### Maestro Protein Preparation Workflow Selected Settings
   1. Convert selenomethionines to methionines
   2. Include Peptides when capping termini.
   3. For HNE, two additional hydrogens on ASN98 and ASN159 had to be added onto the sidechain nitrogen.

1. Cofactor EPE was extracted to its own SDF file cofactor.sdf, other cofactors NAG and FUC were deleted.
2. Using Flare v8. cofactors.sdf was loaded in and was resaved to add appropriate hydrogens not prepared by maestro.
3. Final NMA peptide prepared by maestro was converted to NME by using the builder panel and selecting other edits ... change atom properties.
4. Protein was saved out to protein.pdb

## Ligand Preparation: Steps Taken
1. Utilizing the included edges csv file within the original_structures folder, each ligand found within the included *ligands.sdf file was confirmed to have a node within said edges file. These ligands were then loaded into maestro to confirm their structure and saved to ligands.sdf.
2. Removed the neutral forms of the ligands 16, 19, 21
