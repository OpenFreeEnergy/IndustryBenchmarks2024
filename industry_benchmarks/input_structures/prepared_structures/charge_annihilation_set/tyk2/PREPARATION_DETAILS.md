# Preparation done in Schrodinger suite 2024-2 (patch 1)
# import tyk2_protein.pdb into maestro.  Contains an UNK residue which turns out to be NME.   Remove cap residues
# to allow the prepwizard capping process to work on it. Just used delete atoms, it put them back in the same place anyway.
# Run prepwizard
# Changed residue type of NMA to NME and number to 1178 in the gui
# all done for tyk2
# removed ligands jmc_32neu and ejm_52neu to be able to have at least one charge transformation and no ligand with alternative protonation states.
