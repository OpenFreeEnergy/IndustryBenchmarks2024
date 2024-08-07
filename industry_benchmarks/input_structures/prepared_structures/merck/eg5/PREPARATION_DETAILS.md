Eg5 Target - Merck dataset


Protein:

- In provided input: only standard amino acids, no cofactors, no ions, some water with HOH standard name, termini caps called ACE / NMA

-> No modifications carried out, except NMA -> NME renaming as suggested

- Input validation script was executed with protein file sucessfully


Ligands:

- In provided input: 28 different ligands of which 6 ligands were given in two different protonation states (-NH vs. - NH2+ groups)

- Provided table eg5_extraprotomers_manual_pkacorr_out.csv gave same predicted values for both protonation states

- Published data from Merck did not show which protonation state was used

-> Protonation states with "CHEMBL*_1" - names were removed, leaving 28 different ligands in one protonation state

QUESTIONS:
1. IS THERE A WAY TO DETERMINE WHICH PROTONATION STATES TO USE FOR BENCHMARK?
2. LIGANDS IN SERIES HAVE DIFFERENT TOTAL CHARGES BECAUSE OF -NH+ GROUPS (also in ligand subset that excludes ligands with 2 protonation states), IS THAT INTENTIONAL?


Software used:
Modifications carried out manually directly in files (Vi)


