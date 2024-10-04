
Edits for the jacs_set/thrombin_core system

Notes for adjusting the protein

1.  Use vi to change the resid for 5 residues in chain H:  
ASN 218B -> 150B, 
VAL 247C -> 150C, 
GLY 248D -> 150D, 
LYS 249E -> 150E, 
GLY 250  -> 150F
This is required because the numbering in the original pdb
file made the residue id values nonconsecutive relative to 
the order of the residues along the chain.  This ordering 
is non-pdb compliant and the pdb file would not otherwise
pass the input_validation check.

2. Rename four water molecules with resname SPC to resname HOH.
This was done with Maestro 13.5.128.

3. TYS 63 (chain I) was changed to TYR with removal of the sulfate
group, leaving -OH at the end of the side chain.  This was done with
the 3D builder component of Maestro: delete -SO3 atoms leaving -OH.

4. Cap chain H: GLY 246 (terminating in -CO-H) capped with
NMA 247.  (Use 3D builder to delete terminal H, select C, 
then "Add Fragmane" -> Protein Capping Groups -> NMA) Then
select entire capping group and Optimize just the NMA residue.

5. Cap chain I: TYR 63 (mutated from TYS in step 3 and  terminating 
with -CO-H) capped with NMA 64.  Use 3D builder as above.

6. Cap chain L:  GLU 1C (terminating in -NH2, neutral) was
capped with ACE, as described above.  ILE 14K (terminating in -CO-H)
capped with NMA 15.  Use 3D builder as above.

7. Rename three NMA capping groups to NME using vi.  

Notes for adjsuting the ligands

Ligand file was checked for compounds with multiple conformers and
charge states.  None found; used as supplied.
