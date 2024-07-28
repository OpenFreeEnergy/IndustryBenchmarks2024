  # Preparation done in Schrodinger suite 2024-2 (patch 1)
  # Read in jak1
  # Contains two PTR residues, converted to Tyr with standard Maestro mutate
  # Ran diagnostics, has two orientations for His885.  Select default (53% occ) and commit.
  # this His has the proton on E but is charged on D with no proton.
  # Choose to fix by adding a hydrogen (making it a HIP really).  This maintains charge which I guess is correct
  # Residue name was not modified.
  # The resides from 913 to 916 are missing and the cut-off ends are capped in the input structure
  # Same happens with 946 to 949.   Seem to be missing from multiple structures in the PDB.
  # Obtained fasta from uniprot and added loop and sidechain building to prepwiz settings.   Obtained
  # a structure with modeled loops (protein_loops.pdb).   Suggest running openFE on both systems, capped as provided and with loops.
  # Also had to redo capping of non-loop-fixed version (protein.pdb)
  # All proteins passed validation (SIMULATION COMPLETE message obtained on amazon using CPU to run the input_validation.py script)
  # From the sdf file for ligands, for jak1 I removed the neutral states of ligands with multiple states.
