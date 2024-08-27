
Changes to made to misc/cdk8_koehler system

All edits were made with Maestro 13.5.128, unless noted.

1.  Protein edits

Two chains, A and B:
Chain A N-terminus is LYS 0 ending with -NH2; capped with ACE -1
Chain A C-terminus is GLU 362 ending with -CO-H; capped with NMA 363
Chain B N-terminus is ASP -3 ending with -NH3(+); see below            
Chain B C-terminus is PRO 264 ending with -CO-H; capped with NMA 265
[Capping edits performed using the 3D Builder component of Maestro.  
In each case, select and delete a hydrogen, then select the adjacent atomic
site to the hydrogen, and then "Add Fragment" and select "Protein capping
groups" from the dropdown menus.]

Even though the protocol is to leave charged terminal groups unchanged,
editing had to be done for the N-terminal ASP -3 on the B chain.
The original pdb file had the acid group protonated at the end of the 
ASP side chain (-COOH) and the terminal nitrogen protonated -NH3(+).
This would have called for a NASH residue which is not available in
the Amber 14SB release.  The absence of this residue was preventing
the protein from passing the input_validation check.  Furthermore, 
the experimental conditions that produced the pdb structure were 
nearly neutral (pH 6.9), suggesting the acid group at the end of 
the side chain is possibly deprotonated, -COO(-).  Therefore, the 
proton was removed from the acid group on the ASP side chain that 
was present in the original pdb file, and its adjacent oxygen was 
given a -1 charge, leaving -COO(-).  (Thanks, Irfan)

The  A:ACE capping group was renumbered from 1 to -1 using the 
"Other Edits" and "Change Atom Properties" component of 3D Builder.
The edited structure was exported with "Reorder by residue number"
option checked, which properly places the capping groups at the 
start/end of the chains.

vi was used to change "NMA" to "NME" in the pdb file.

Alternative conformations existed for 6 residues, highest occupancy was 
selected:  A:CYS 64, A:ARG 178, A:HIS 248, B:ARG 107, B:ARG 247, 
B HOH 438.  Note that for B:ARG 247 the occupancy of both conformers 
was exactly 50% and the first (A) conformer was selected.

Other than noted above for B:ASP -3, protonation states were checked 
and seen to be unchanged. No unaccounted changes in charge sites or 
number of atoms and residues.


2.  Ligand edits

Ten ligands were checked to see that there was only one conformer per 
compound and only one protonation/charge state per compound.  All compounds
were neutral (uncharged).

