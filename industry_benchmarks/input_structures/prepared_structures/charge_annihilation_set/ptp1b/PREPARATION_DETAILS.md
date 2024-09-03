## System preparation details

* Original files were downloaded from:
	- Protein: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/main/industry_benchmarks/input_structures/original_structures/charge_annhil/ptp1b_protein.pdb
	- Ligands: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/main/industry_benchmarks/input_structures/original_structures/charge_annhil/ptp1b_ligands.sdf

# Protein preparation
* The atom name "O   ASP A 298" was manually set to "O2   ASP A 298" (line 4808 in protein.pdb)
* The protein.pdb passed the validation (python input_validation.py --pdb protein.pdb)

# Ligand preparation
* The ligand 26neu was removed
