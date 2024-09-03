## System preparation details

* Original files were downloaded from:
	- Protein: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/main/industry_benchmarks/input_structures/original_structures/charge_annhil/ephx2_protein.pdb
	- Ligands: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/main/industry_benchmarks/input_structures/original_structures/charge_annhil/ephx2_ligands.sdf

# Protein preparation
* ACE and NMA caps were added using the "build" function in Schrodinger's Maestro (Release 2023-4)
	- Final residue numbers: ACE236 and NMA549
* NMA was manually renamed to NME
* Hydrogen names of ACE and NME were manually altered 
* The protein.pdb passed the validation script (python input_validation.py --pdb protein.pdb)

# Ligand preparation
* Ligands 47neu and 19neu were removed
