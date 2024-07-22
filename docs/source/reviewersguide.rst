*********************************
Guide for reviewers / maintainers
*********************************

This page outlines key instructions for those who are reviewing and/or
maintaining various parts of the industry benchmark process. Under most
circumstances, this responsibility should be limited to the OpenFE staff
members.

Please see individual headings for instructions for a specific task.


Public Datasets Phase 1: Evaluating System Submissions
======================================================

All submitted files should follow the instructions provided under the
:ref:`contributing inputs <contributing-inputs>` page. The contributors
will have a preparation checklist for each contribution. The main aim of
reviewing is ensuring that this checklist has been appropriately completed.

Below are a series of things you should consider when reviewing new
pull requests:


1. How the data is organized
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **directory structure** and **file names** should match those given
in the `templates directory <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/reviewers-guide/industry_benchmarks/input_structures/prepared_structures/template>`_.


2. Contents of `PREPARATION_DETAILS.md`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A `PREPARATION_DETAILS.md` file should have been created for each system.
This file should outline how the the systems were prepared.

**Note that this is what we will use to generate the methods
information of the paper** . There should be sufficient detail so
that you could reproduce the work done should it be necessary.

You should consider:
  - Do the contents of the file make sense?
  - Are there version numbers for any tools employed?
  - Is there a need to justfify any decisions (i.e. if a choice was made to
    change the inputs from those provided by Schrodinger, in a different
    manner to what we ask in the preparation instructions, this needs to
    be explained).
  - Do the details provided differ from what was done to the files?


3. Checking the ligands
~~~~~~~~~~~~~~~~~~~~~~~

A `ligands.sdf` file should have been created for each system which contains
**only one** ligand state.

To verify this, you can use the `ligand prep checking script <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/reviewers-guide/industry_benchmarks/utils/maint/check_ligand_prep.py>`_.

Please note that this script is not infallible, you should have a quick look
through the file for the following:
  - Are there any entries with the same ligand name but additional modifiers (e.g. neu or chg or R)?
  - If you load the file into pymol with the protein, do the ligands sit in the binding site?


4. Checking the protein
~~~~~~~~~~~~~~~~~~~~~~~

A `protein.pdb` file should have been created for each system.

You should check that the PDB has been prepped according to the instructions.
To help you with this, you can use the `protein prep checking script <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/reviewers-guide/industry_benchmarks/utils/maint/check_protein_prep.py>_`.

This script will help you do the following:
  - Check that the number of atoms in each residue (excluding the termini)
    does not change (i.e. protonation states have not changed).
  - Check that the positions of protein atoms have not change significantly.
  - Check that the number of non-protein atoms match (note this could be ok
    if cofactors were removed from the file or caps were renamed from UNK).
  - Check for any disulfide bridges and, if there are, make sure that the
    cysteine SG atoms do not have a hydrogen bound (ideally with a CONECT record for the SG-SG bond).

You should also manually check the following:
  - That the termini are capped / uncapped based on how Schrodinger handled the original PDB files.
    If the termini are in different capping states then the contributor should have a scientific reason why.
  - Any cofactors were removed to a separate `cofactors.sdf` file and has appropriate bond order / protonation states.
  - Crystallographic waters and metals were kept.
  - Any glycosylation sites have been reverted to normal AAs (with the sugars removed).
  - Any PTMs have been reverted to their canonical AA counterparts.


5. Ensuring the MD validation script runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whilst each contributor should have run the
:ref:`validation script <input-validation>`, it is worth double checking that
things still run, so please do run the script on all files provided before
merging.

