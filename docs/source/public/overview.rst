.. _public_overview:

*************************
Public Dataset Benchmarks
*************************


This page outlines the plans and the instructions for the public dataset portion of the 2024 OpenFE Industry Benchmark study.


Overview
********


The Public Benmarking Study concentrates on validating the use of the `OpenFE Toolkit <https://docs.openfree.energy/en/stable/>`_
on publicly available protein-ligand binding datasets. Specifically we concentrate
on re-calculating the `Schrodinger public binding free energy benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark>`_
from the `2023 large scale study by Ross et al. <https://www.nature.com/articles/s42004-023-01019-9>`_.


.. figure:: public_overview.png
   :scale: 100%
   :align: center
   :alt: Figure 1 - Overview of Public Dataset Benchmark process.


*Figure 1: Overview of Public Dataset Benchmark process.*

As outlined in Figure 1, the public benchmark will be split into four phases:

1. Remediation and validation of input structures from the `Schrodinger public binding free energy benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark>`_.
2. Simulating relevant free energy perturbations on each participating industry member’s own HPC resources.
3. Extracting and gathering results.
4. Processing results into a manuscript for submission to an appropriate journal.

At each stage of the process, the OpenFE team will provide necessary inputs and support to all participating members.


.. _public_phase1:

Phase 1: Preparing Inputs
*************************

In this phase, benchmark input structures will be prepared by industry partners for use with the OpenFE toolkit.


**Start date:** *Early May 2024*

**End date:** *Early June 2024*


Allocation Benchmark Systems
============================

Each participating industry partner will be expected to select one or more system for benchmarking.

Please see the :ref:`benchmark system overview <benchmark_systems>` for a list of all systems and which
of these have already been allocated.

At a very minimum, we aim to simulate these entirety of the following three datasets:
* JACS
* Fragments
* Merck

Please :ref:`get in touch <get_in_touch>` with the OpenFE team should you wish to select
any new systems or if you wish to be reminded of which systems you have been allocated.


Input Data Source
=================

Input structures, perturbation networks, and FEP+ results are be taken from the `v2.0 release of the Schrodinger 2023 Public Benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/v2.0>`_.

For convenience, a `snapshot of these inputs <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/input_structures/original_structures>`_ have been provided in the OpenFE 2024 benchmark repository.


Remediation of inputs
=====================

Input structures will need to be adapted in order to be used with the OpenFE toolkit.

Nominally input remediation will involve:

* Extracting any cofactors from input PDB files
* Fixing capping groups on terminal residues
* Mutating any non-canonical amino acids
* Addressing any issues with residue and atom names

Depending on the software used to prepare the inputs, the resulting files may need to be stripped of extraneous or sensitive metadata.

A list of instructions and conditions for preparing inputs are provided :ref:`here <input_preparation>`.
Updates will be made to these instructions based on feedback by benchmark partners.


.. toctree::
   :maxdepth: 1

   input_preparation


Depositing inputs
=================

All remediated inputs will be deposited in the `OpenFE 2024 benchmark repository <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024>`_.
This will allow the OpenFE team to:

1. Gather preparation conditions to be included in relevant publications
2. Check-in with industry partners and gather feedback on the the input preparation experience
3. Allow for the OpenFE team to help with any unanticipated issues


Please see the :ref:`contributing inputs instructions <contributing-inputs>` for more details on this process.


.. toctree::
   :maxdepth: 1

   contributing_inputs


Phase 2: Running Simulations
****************************


.. note::
   Details for phase 2 of the public dataset benchmarks are still being finalised. These will be updated as soon as possible!


In this phase, industry partners will run alchemical transformations for their allocated systems on their HPC resources.

**Start date:** *Early June 2024*

**End date:** *End of August 2024*

Please note that we expect the private dataset industry benchmark to start alongside this phase.


Simulation Planning: LOMAP networks
===================================

*Details to be published very soon!*


Simulation execution
====================

All planned simulations will be run by industry partners on their own clusters using OpenFE execution tooling,
i.e. through the `quickrun method <https://docs.openfree.energy/en/latest/guide/execution/quickrun_execution.html>`_.


Compute Requirements
====================

The following compute resources will be required:

**GPU Hardware**

Industry partners are expected to have the following GPU hardware:

* Approximately 24 GPU hours per triplicate repeat of each standard transformation
   * Up to 15 GPU days for net charge transformations
* CUDA 10.2 or above
* Non-exclusive compute mode
* Assignment of a single GPU ID per openfe quickrun execution (i.e. by setting CUDA_VISIBLE_DEVICEs if necessary)

**Data storage**

Industry partners will be expected to keep simulation outputs for the duration of the study, in case the data needs to be post-processed during the publication stage.

We estimate a requirement of **5 GB per alchemical transformation** edge.


Phase 3: Data Analysis
**********************


.. note::
   Details for phase 3 of the public dataset benchmarks are still being finalised. These will be updated as soon as possible!


In this phase, relevant simulation results will be gathered from industry partners.

**Start date:** *Early September 2024*

**End date:** *End of October 2024*


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

Industry partners will be expected to submit this information back to OpenFE for
analysis. Please note that all data will be collected in a human readable format
in order to allow industry partners the ability to review the data ahead of submission
back to the OpenFE team.

Analysis of results
===================

**Analysis of individual systems**

Initial analysis of results for each system will be carried out by each industry benchmark partner
with the help of the OpenFE development team.

**Analysis of all results**

A final analysis of all simulation results will be conducted by the OpenFE development team with
help from volunteering industry board and technical advisory committee members.


Phase 4: Paper writing
**********************


.. note::
   Details for phase 3 of the public dataset benchmarks are still being finalised. These will be updated as soon as possible!


**Start date:** *Early November 2024*

**End date:** *December 2024*


Drafting the manuscript
=======================

Once all results have been gathered, the OpenFE team alongside volunteer members of
the OpenFE board and technical advisory committee will draft a manuscript for open access
publication at a relevant journal.

Review of manuscript and authorship
===================================

All authors will be expected to review and approve the manuscript prior to journal submission.
We anticipate doing this in a two round process, the first round where authors are invited
to comment on the manuscript, followed by a second one for legal review by each partner organization.

Authorship will be offered to all those involved in the benchmarking process and inclusion will be
left to the discretion of each organization.

