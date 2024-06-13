.. _input-preparation:

*****************************************************************
Instructions for preparing inputs for the OpenFE public benchmark
*****************************************************************

Overview & Aims
***************

This page provides guidance on preparing the `Schrodinger benchmark files <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/v2.0/fep_benchmark_inputs/structure_inputs>`_
for use in the OpenFE 2024 public dataset industry benchmark.

Whilst each partner organization will be responsible for preparing inputs
using their own tooling, this guidance aims to provide a set of rules to
standardize the process and ensure that scientifically equivalent simulation inputs are generated.

This page will be continually updated based on feedback from benchmark partners.
Please :ref:`reach out <get_in_touch>` should you have any questions or
require additional information whilst preparing your inputs.


Checklist
*********

By the end of these instructions you should have:

* A PDB file containing the system protein, crystallographic waters and ions (named `protein.pdb`)
* An SDF file containing the ligands being mutated (named `ligands.sdf`)
* (Optional) An SDF file containing any system cofactors (named `cofactors.sdf`)


Input preparation instructions
******************************

Input data source
=================

All input data are sourced from the `v2.0 release of the Schrodinger binding free energy benchmarks <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/v2.0>`_.

For your convenience, a snapshot of these inputs have been placed under the
OpenFE Public Benchmark repository. Specifically, input PDB, SDF and edges
CSV files can be found under the relevant `original structures <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/structure_inputs/original_structures>`_
sub-folder.

We recommend that you clone this repository and use those input files.

Extracting cofactors
====================

In its benchmark inputs, Schrodinger places cofactor molecules within the PDB file.
Currently OpenFE cannot parse such inputs as small molecules are expected as separate inputs.

Should your input PDB have cofactors, please extract these from the PDB and store them
within an SDF file named ``cofactors.sdf``. You may be required to manually assign and/or
correct bond orders such that they represent the expected state of the molecules and the
SDF files are readable by RDKit.

.. note::
   Waters and ions do not need to be extracted from the PDB file.

Capping proteins
================

Whether or not proteins have to be capped depends on how the protein’s termini are handled
in the Schrodinger provided PDB files.

1. If the termini are neutral, e.g. NH2 or C(=O)H, or have ACE/NME caps:
    * Assign ACE and/or NME caps to the termini
.. note::
   In some cases Schrodinger / Maestro may call the NME cap NMA, if this happen, the cap should be renamed to NME.

2. If the termini are charged, e.g. COO- or NH3+:
    * Keep the termini in their charged state
.. note::
   In cases where this is observed, e.g. JNK1, there is a possibility that interactions between the perturbed ligand and the termini may form.
   Keeping the termini charged should retain the intended interaction.

.. note::
   Some tools, e.g. Pymol, are known to erroneously place caps out of order with the protein chain or TER cards between the protein chain and the cap. Please ensure that the caps appear in order (i.e. before or after the relevant termini residue), and that there are no TER cards between the termini residue and the cap.

Non canonical amino acids (PTMs)
================================

Some of the input structures have non-canonical amino acids, specifically TPO, PTR, or TYS residues.
The OpenFE team has reviewed all such cases and the amino acids were found to be far from the binding site.

As the OpenFE software cannot easily handle PTMs at this stage, **these residues should be modified back to their canonical alternative**.
Tools such as Pymol offer the ability to mutate residues in this manner.

Protonating inputs
==================

All files in the Schrodinger input files are considered to be in their intended protonation
state. No additional protonation steps should be carried out.

.. note::
   1. Care should be taken that any processing steps do not alter the protonation state. For example, PyMol is known to change histidines from HID to HIE during capping.

   2. To ensure that protonation states were not altered, please ensure that the prepared PDB file has the same number of protein atoms outside of capping groups.

Setting residue names
=====================

Where possible, residues should be assigned PDB-compliant names.

*Example 1: Waters named SPC (e.g. in the case of Thrombin in the JACS set), should be renamed to HOH.*

*Example 2: Capping groups named NMA should be renamed to NME (e.g. in the case of PTP1B in the JACS set).*

Fixing hydrogen atom names
==========================

In some cases, hydrogen names may need to be manually altered to match expected, i.e. PDB compliant, names.

These exact cases can be difficult to identify, running the validation script (see below), will help identify these. Please reach out to the OpenFE team should you encounter any unknown hydrogen names.

*Example 1: GLY termini hydrogens being named 3HA and HA instead of HA3 and HA2.*

*Example 2: HIS (in the HID state) hydrogens being named 1HD, 2HD, and 1HE instead of HD1, HD2, and HE1.*

Validating prepared files
=========================

To ensure that prepared files can be run using OpenFE, a short MD simulation validation script has been provided under
`utils/input_validation.py <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/utils/input_validation.py>`_.
In an environment with OpenFE 1.0 installed, please run this script by calling:

.. code-block:: python

   # If you don’t have cofactors
   python input_validation.py --pdb protein.pdb

   # If you have cofactors
   python input_validation.py --pdb protein.pdb --cofactors cofactors.sdf


If the script outputs “SIMULATION COMPLETE”, then your inputs are suitable for use with OpenFE. If they do not, then there is likely an issue with the input file. Please report the error message emitted when contacting the OpenFE team for advice on how to fix any issues.

.. note::
   This script runs a very short simulation, it is recommended that it is executed on a machine with a CUDA-enabled GPU.

Preparing the ligand file
=========================

For some datasets, the Schrodinger public binding free energy benchmark set includes multiple binding modes (e.g. different rotamers) 
and protonation states of ligands. For this current study, we will only consider a single conformation and protonation state for each of the ligands. 

If the dataset contains ligands in multiple conformations or protonation states, the state that likely contributes the most to binding should be identified (by looking at previous results) and the less favorable state should be removed from the input ``ligands.sdf`` file.

The FEP+ ligand predictions can be found `here <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/21_4_results/ligand_predictions>`_.

In the following, we will go into the details on how to extract the necessary information for ligands with multiple binding modes, multiple protonation states, and multiple stereo isomers.

**1. Multiple binding modes**

For ligands that were run in multiple binding modes, the table of FEP+ ligand predictions reports only the binding mode
that was calculated to contribute more to binding.

*Example: JNK1 (JACS set)*

* Opening the `Table of ligand predictions <https://github.com/schrodinger/public_binding_free_energy_benchmark/blob/main/21_4_results/ligand_predictions/jacs_set/jnk1_manual_flips_symbmcorr_out.csv>`_
* The table shows the experimental and calculated binding free energies for 21 ligands, while there had been 38 nodes in the FEP+ network
* Remove all ligands from the ``ligands.sdf`` file that are not listed in this table
* e.g. ``18637-1`` is present in the table but not ``18637-1 flip``, therefore we would remove ``18637-1 flip``
* It may also be helpful to look at the `Table of edge predictions <https://github.com/schrodinger/public_binding_free_energy_benchmark/blob/main/21_4_results/edge_predictions/jacs_set/jnk1_manual_flips_out.csv>`_
  to identify the ligand pairs for which multiple binding modes had been used
* e.g. first edge between ligand ``18637-1`` and its alternate binding mode ``18637-1 flip``

**2. Multiple protonation states**

For ligands for which multiple protonation states were included in the ligand network,
the table of FEP+ ligand predictions reports calculated binding free energies from all states.
The values include a pka correction as described in work by `Oliveira et al <https://pubs.acs.org/doi/10.1021/acs.jctc.8b00826>`_.
For this study we will be using a single protonation state per ligand, choosing the protonation state that had been used in the original studies by `Schindler et al. (Merck set) <https://pubs.acs.org/doi/10.1021/acs.jcim.0c00900>`_, `Chen et al. (charge annihilation set) <https://pubs.acs.org/doi/10.1021/acs.jctc.8b00825>`_, and `Cappel et al. (MCS docking set) <https://pubs.acs.org/doi/10.1021/acs.jcim.9b01118>`_.

From the systems picked by industry partners, as of writing these instructions, the following systems have ligands in multiple protonation states:

* Merck set: EG5, TNKS2
* MCS docking set: HNE
* Charge annihilation set: JNK1, EGFR, DLK, JAK1, TYK2, ITK, CDK2

If you picked one of these systems, please reach out to us with any questions regarding the protonation state assignment!
 

**3. Multiple stereo isomers**

For ligands where multiple stereo isomers where included in the ligand network,
the table of FEP+ ligand predictions reports results from both stereo isomers.
In this case we will keep the stereo isomer with the more negative calculated binding free energy and remove the other stereo isomer from the ``ligands.sdf`` file.

*Example: MUP-1 (fragments dataset)*

* Opening the `Table of ligands predictions <https://github.com/schrodinger/public_binding_free_energy_benchmark/blob/main/21_4_results/ligand_predictions/fragments/frag_mup1_out.csv>`_
* For ligand ``SBT`` there are results from two stereo isomers, ``SBT_R`` and ``SBT_S``
* The calculated binding free energy of ligand ``SBT_S`` is more negative than for ligand ``SBT_R`` (-9.14 vs. -8.77 kcal/mol)
* In this case we would remove ligand ``SBT_R`` from the ``ligands.sdf`` file


Submitting prepared input files
===============================

All prepared inputs should be submitted to the OpenFE Public Benchmark github repository, more specifically to the
`prepared_structures <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/inputs/prepared_structures>`_ subfolder.
This should be done via Pull Request, with a folder for each prepared system including the protein PDB, ligand SDF, and if available cofactor SDF file.
A short bullet point summary of any remediation steps, including any software used, should also be included as a markdown file.
Further details can be found in the :ref:`contributing-inputs` page.

If necessary, you may email the OpenFE team with this information and the Pull Request will be opened on your behalf.

Once the Pull Request is opened, the OpenFE team will carry out a minimal review of the contents, including a short validation that the alchemical transformations will work. If all checks pass, the Pull Request will be merged and you should be ready to start the next step in the benchmarking process (setting up the alchemical network).
