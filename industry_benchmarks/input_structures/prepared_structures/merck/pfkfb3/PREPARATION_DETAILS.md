Prepared by Hendrik Geoddeke

- Import original inputs from the repo in Maestro 2024-1
- Check the sequence to see if it contains any non-natural amino acids
- Check all caps to see if any need to be charged or not.
- Open the Protein preparation wizard (Schrodinger 2024-1)and run only step 3 - "Preprocess". The only three options activated are:
	->Cap termini
	->Convert selomethionines to methionines
	->Include peptides when capping termini
- Change the cap names if not NME by using Schrödingers "Builder"
- In case of multiple rotameric states for the same ligand, the one with the most negative dG is selected and the others are discarded.
- Export outputs, named "ligand.sdf", "protein.pdb" and "cofactors.sdf"
- Run input through validation script

*For some reason, all three cofactors had a strange Lewis structure, where all oxygens with a double bond were negatively charged(e.g., C=O-). The Lewis structures seem to be correct in the FEP+ output (pfkfb3_automap_symbmcorr_out.fmp) and the states of the cofactors were changed manually in Maestro according to the pfkfb3_automap_symbmcorr_out.fmp file.
