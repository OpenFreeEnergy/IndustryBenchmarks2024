
Changes made to merck/cdk8_5cei_new_helix_loop_extra.

Unless noted, all edits were performed  with Maestro 13.5.128, unless noted.

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
The original pdb file had the acid group at the end of the ASP side chain
protonated (-COOH) and the terminal nitrogen protonated -NH3(+).
This would have called for a NASH residue which is not available in
the Amber 14SB release.  The absence of this residue was preventing
the protein from passing the input_validation check.  Furthermore, 
the experimental conditions that produced the pdb structure were 
nearly neutral (pH 6.9), suggesting the acid group at the end of 
the ASP side chain is possibly deprotonated, -COO(-).  Therefore, the 
proton was removed from the acid group on the ASP side chain that 
was present in the original pdb file, and its adjacent oxygen was 
given a -1 charge, leaving -COO(-). This edit allowed the protein
to pass the input_validation check.   (Thanks, Irfan)

The  A:ACE capping group was renumbered from 1 to -1 using the 
"Other Edits" and "Change Atom Properties" component of 3D Builder.
The edited structure was exported with "Reorder by residue number"
option checked, which properly places the capping groups at the 
start/end of the chains.

vi was used to change "NMA" to "NME" in the pdb file.

Other than noted above for B:ASP -3, protonation states were checked 
and seen to be unchanged. No unaccounted changes in charge sites or 
number of atoms and residues.


2.  Ligand edits

The original sdf file had 41 ligand entries, consisting of 33 distinct
compounds, 8 of which had two docking poses.  Only one each of these
eight pairs of docking poses should be included.  These are the poses included
in the ligand_predictions/merck/cdk8_5cei_new_helix_loop_extra directory:
13 13-flip:  keep 13
16 16-flip:  keep 16-flip
17 17-flip:  keep 17-flip
18 18-flip:  keep 18-flip
42 42-flip:  keep 42-flip
43 43-flip:  keep 43
44 44-flip:  keep 44-flip
45 45-flip:  keep 45-flip
After these deletions, there woud be  33 participating ligands, however,
ligand 28 was not included in the Schrodinger ligand_prediction file,
and this was also omitted.  The final count is 32 participating ligands:
13      18-flip 23     29     34     39      44-flip
14      19      24     30     35     40      45-flip
15      20      25     31     36     41
16-flip 21      26     32     37     42-flip
17-flip 22      27     33     38     43
 
Changed the names of the six "flipped" ligands in the ligands.sdf file 
to add a hyphen instead of a space so that file names generated from 
the ligand names by OpenFE software will not have a space. 
E.g., "13 flipped" was changed to "13-flipped".

The ligands were checked to see that there was only one conformer per 
compound and only one protonation/charge state per compound.  All compounds
were neutral (uncharged).

