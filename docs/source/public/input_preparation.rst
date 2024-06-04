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
* A CSV file containing the original FEP+ perturbation network (named `edges.csv`)
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

For the ligand inputs, the Schrodinger public binding free energy benchmark set includes multiple binding modes (e.g. different rotamers) 
and protonation states for some of the ligands. For this current study, we will only consider a single conformation and protonation state of the ligands. 
Here, we will be using the binding mode and protonation state that was predicted to be most potent by FEP+.

If the dataset contains ligands in multiple conformations or protonation states, the more favorable state should be identified and the less favorable state removed from the input ``ligands.sdf`` file.

**Assessing the more potent binding mode or protonation state** 

The FEP+ edge predictions can be found `here <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/21_4_results/edge_predictions>`_. 

Example: JNK1 (JACS set)

* `Table of edge predictions <https://github.com/schrodinger/public_binding_free_energy_benchmark/blob/main/21_4_results/edge_predictions/jacs_set/jnk1_manual_flips_out.csv>`_
* For all edges that connect different binding modes of the same ligand, the experimental ddG values have a value of 0.0 kcal/mol. 
* The calculated ddG value between the first edge (ligand *18637-1* and its alternate binding mode *18637-1 flip*) is 2.54 kcal/mol. 
* This means that the original binding mode (*18637-1*) is predicted to be more favorable than the flipped binding mode (*18637-1 flip*). 
* In this case we would remove ligand *18637-1 flip* from the dataset.

Submitting prepared input files
===============================

All prepared inputs should be submitted to the OpenFE Public Benchmark github repository, more specifically to the
`prepared_structures <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/inputs/prepared_structures>`_ subfolder.
This should be done via Pull Request, with a folder for each prepared system including the protein PDB, ligand SDF, relevant edges CSV, and if available cofactor SDF file.
A short bullet point summary of any remediation steps, including any software used, should also be included as a markdown file.
Further details can be found in the :ref:`contributing-inputs` page.

If necessary, you may email the OpenFE team with this information and the Pull Request will be opened on your behalf.

Once the Pull Request is opened, the OpenFE team will carry out a minimal review of the contents, including a short validation that the alchemical transformations will work. If all checks pass, the Pull Request will be merged and you should be ready to start the next step in the benchmarking process (setting up the alchemical network).
