.. _phase_1_4_instruction:

Instructions for the OpenFE 2024 Public Benchmark
#################################################

This page outlines the plans and the instructions for the four phases of the public dataset portion of the 2024 OpenFE Industry Benchmark study.
The work outlined here will be undertaken by either participating industry members or the OpenFE team.

Overview
********

TODO: Add Figure from Irfan

As outlined in Figure 1, the public benchmark will be split into four phases:

1. Remediation and validation of input structures from the `Schrodinger public binding free energy benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark>`_.
2. Simulating relevant free energy perturbations on each participating industry member’s own HPC resources.
3. Extracting and gathering results.
4. Processing results into a manuscript for submission to an appropriate journal.

At each stage of the process, the OpenFE team will provide necessary inputs and support to all participating members.

Phase 1: Benchmark inputs
*************************

In this phase, benchmark input structures will be prepared by industry partners.

Timeline
========

**Expected start date.**
Industry partners are expected to start work in preparing inputs by the **end of April 2024**.

**Expected completion date.**
We expect partners to complete this phase by the **end of May 2024**.
Please note that industry partners do not need to wait on others to finish this task before they carry on with phase 2.

Input Data Source
=================

Input structures, perturbation networks, and FEP+ results will be taken from the `v2.0 release of the Schrodinger 2023 Public Benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/v2.0>`_.

This dataset accompanies the 2023 manuscript `“The maximal and current accuracy of rigorous protein-ligand binding free energy calculations” <https://www.nature.com/articles/s42004-023-01019-9>`_ by Ross et al.

Remediation of inputs
=====================

Whilst Schrodinger’s team has kindly provided all simulation inputs to their benchmarks, they unfortunately cannot be directly used with various different open source tools, including those employed by OpenFE. To resolve this, industry partners will need to remediate these inputs as necessary.

Nominally input remediation will involve:

* Extracting any cofactors from input PDB files
* Fixing capping groups on terminal residues
* Mutating any non-canonical amino acids
* Addressing any issues with residue and atom names

Depending on the software used to prepare the inputs, the resulting files may need to be stripped of extraneous or sensitive metadata.

A list of instructions and conditions for preparing inputs are provided :ref:`here <input_preparation>`.
Updates will be made to these instructions based on feedback by benchmark partners.

Depositing inputs
=================

All remediated inputs will be deposited in the `OpenFE 2024 benchmark repository <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024>`_.
This will allow the OpenFE team to:

1. Gather preparation conditions to be included in relevant publications
2. Check-in with industry partners and gather feedback on the the input preparation experience
3. Allow for the OpenFE team to help with any unanticipated issues

Phase 2: Running Simulations
****************************

In this phase, industry partners will run alchemical transformations on their HPC resources.

Timeline
========

**Expected start date.**
Industry partners are expected to start simulations by the **end of May 2024**.

**Expected completion date.**
We expect partners to complete this phase by the **end of August 2024**.
Please note that we expect the private dataset industry benchmark to start alongside this phase.

Simulation Planning: LOMAP networks
===================================

To be added after edited in the planning document

Simulation execution
====================

All planned simulations will be run by industry partners on their own clusters using OpenFE execution tooling, i.e. through the `quickrun method <https://docs.openfree.energy/en/latest/guide/execution/quickrun_execution.html>`_.

**Expected compute requirements**

The following compute resources will be required:
**GPU Hardware**

Industry partners are expected to have the following GPU hardware:

* Approximately 24 GPU hours per triplicate repeat of each standard transformation
   * Up to ??? GPU days for net charge transformations
* CUDA 10.2 or above
* Non-exclusive compute mode
* Assignment of a single GPU ID per openfe quickrun execution (i.e. by setting CUDA_VISIBLE_DEVICEs if necessary)

**Data storage**

Industry partners will be expected to keep simulation outputs for the duration of the study, in case the data needs to be post-processed during the publication stage.

We estimate a requirement of **5 GB per alchemical transformation** edge.

Phase 3: Results Analysis and Gathering
***************************************

In this phase, relevant simulation results will be gathered from industry partners.

Timeline
========

**Expected start date.**
Gathering of simulation results is expected to begin as soon as possible, but no later than the start of **September 2024**.

**Expected completion date.**
We expect partners to complete this phase by the end of **October 2024**.

Gathering of results
====================

Industry partners will be expected to post-process simulation outputs using a specialized script provided by the OpenFE team.

This script will:

* Extract relevant free energy estimates, including time series of free energies
* Gather simulation health metrics
   * Overlap matrix and replica exchange probability plots
   * Relevant structural analysis plots
* Gather additional simulation information (optional)
   * Additional simulation metrics, relevant for the OpenFE 2024 scoring data project, may be gathered.

Industry partners will be expected to submit this information back to OpenFE for analysis. Please note that all data will be collected in a human readable format in order to allow industry partners the ability to review the data ahead of submission back to the OpenFE team.

**Estimated development cost**

Development of the necessary analysis scripts and their documentation is expected to take **2 FTE weeks**.

Analysis of results
===================

**Analysis of individual systems**

Initial analysis of results for each system will be carried out by each industry benchmark partner with the help of the OpenFE development team. Should any issues be identified, additional work in data gathering and/or simulations may be required.

**Analysis of all results**

A final analysis of all simulation results will be conducted by the OpenFE development team with help from volunteering industry board and technical advisory committee members.

**Time investment**

We estimate this task to require an estimated **2 FTE months** of OpenFE developer time.

This includes:

1. The development of specialized scripts to analyze perturbation networks containing multiple copies of the same ligand in different conformational and protonation states.
2. Time spent with industry partners investigating non-ideal simulation results
3. Time spent gathering results and creating appropriate meta analyses and plots

Phase 4: Paper writing
**********************

Timeline
========

**Expected start date.**
Drafting of the manuscript is expected to start in **November 2024**.

**Expected completion date.**
We expect partners to review the manuscript in **December 2024**.

Drafting the manuscript
=======================

Once all results have been gathered, the OpenFE team alongside volunteer members of the OpenFE board and technical advisory committee will draft a manuscript for open access publication at a relevant journal.

Review of manuscript and authorship
===================================

All authors will be expected to review and approve the manuscript prior to journal submission. We anticipate doing this in a two round process, the first round where authors are invited to comment on the manuscript, followed by a second one for legal review by each partner organization.

Authorship will be offered to all those involved in the benchmarking process and inclusion will be left to the discretion of each organization.








