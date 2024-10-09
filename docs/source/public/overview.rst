.. _public_overview:

*************************
Public Dataset Benchmarks
*************************


This page outlines the plans and the instructions for the public dataset portion of the 2024 OpenFE Industry Benchmark study.


Overview
********


The Public Benchmarking Study concentrates on validating the use of the `OpenFE Toolkit <https://docs.openfree.energy/en/stable/>`_
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


Before you start
****************

.. warning::
   Please do not run benchmarks with pre-existing installations of OpenFE.
   You should look at our specific :ref:`benchmark openfe install instructions <installation>`
   for how to set up openfe.


Before starting you should :ref:`install openfe <installation>` on any local machine you will use
to prepare the benchmark systems and any HPC systems you will be running simulations on.

Please note that a very specific install of `openfe` and its dependencies is required in
order to reduce any variance introduced by new software releases throughout the length
of the study.


.. _public_phase1:

Phase 1: Preparing Inputs
*************************

In this phase, benchmark input structures will be prepared by industry partners for use with the OpenFE toolkit.


**Start date:** *Mid June 2024*

**End date:** *Mid July 2024*


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

Input structures, perturbation networks, and FEP+ results are taken from the `v2.0 release of the Schrodinger 2023 Public Benchmark set <https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/v2.0>`_.

For convenience, a `snapshot of these inputs <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/input_structures/original_structures>`_ has been provided in the OpenFE 2024 benchmark repository.


Remediation of inputs
=====================

Input structures will need to be adapted in order to be used with the OpenFE toolkit.

Nominally input remediation will involve:

* Extracting any cofactors from input PDB files
* Fixing capping groups on terminal residues
* Mutating any non-canonical amino acids
* Addressing any issues with residue and atom names
* Removing alternate ligand conformational or protonation states

Depending on the software used to prepare the inputs, the resulting files may need to be stripped of extraneous or sensitive metadata.

A list of instructions and conditions for preparing inputs are provided :ref:`here <input-preparation>`.
Updates will be made to these instructions based on feedback by benchmark partners.


.. toctree::
   :maxdepth: 1
   :hidden:

   input_preparation


Depositing inputs
=================

All remediated inputs will be deposited in the `OpenFE 2024 benchmark repository <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024>`_.
Doing this will involve a lightweight review process by the OpenFE team, which will allow them to:

1. Gather preparation conditions to be included in relevant publications
2. Check-in with industry partners and gather feedback on the input preparation experience
3. Allow for the OpenFE team to help with any unanticipated issues


Please see the :ref:`contributing inputs instructions <contributing-inputs>` for more details on this process.


.. toctree::
   :maxdepth: 1
   :hidden:

   contributing_inputs


.. _public_phase2:

Phase 2: Running Simulations
****************************


.. warning::
   Please do not start Phase 2 until your prepared benchmark inputs have been approved by the OpenFE team.


In this phase, industry partners will run alchemical transformations for their allocated systems on their HPC resources.

**Start date:** *Early July 2024*

**End date:** *End of August 2024*

Please note that we expect the private dataset industry benchmark to start alongside this phase.

.. _simulation_planning:

Simulation Planning: LOMAP networks
===================================

The setup of the relative binding free energy network is carried out using the planning script provided under 
`utils/plan_rbfe_network.py <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/utils/plan_rbfe_network.py>`_.

This script will carry out the following steps:

* Loading the ligands, protein, and (if provided) cofactors
* Computing partial charges for ligands and cofactors using antechamber AM1BCC
* Creating a network of ligand transformations using the Kartograf atom mapper, LOMAP scorer, and LOMAP network generator
* Assigning settings to the transformations. Settings differ depending on whether the ligand transformation involves a change in net charge
   * Non charge changing transformations: 11 lambda windows, 5 ns production run per lambda window
   * Charge changing transformations: 22 lambda windows, 20 ns production run per lambda window
* Creating the ``AlchemicalTransformation``\ s for solvent and complex legs and saving them to disc as json files

In an environment with OpenFE 1.0 installed, please run this script by calling:

.. code-block:: bash

   # If you don’t have cofactors
   python plan_rbfe_network.py --pdb protein.pdb --ligands ligands.sdf --output network_setup

   # If you have cofactors
   python plan_rbfe_network.py --pdb protein.pdb --ligands ligands.sdf --cofactors cofactors.sdf --output network_setup

This command will create a folder (named ``network_setup`` as specified using the ``--output`` flag) that contains a separate ``.json`` file for the solvent and complex legs 
for every edge in the network. The folder also contains a ``ligand_network.graphml`` file that is a serialized version of the ``LigandNetwork``.

.. warning::
   Since the partial charge assignment can be slow, we recommend putting the planning command in a bash script and executing it on a high performance workstation or HPC resource. 

.. _simulation_execution:

Simulation Execution
====================

All planned simulations will be run by industry partners on their own clusters using OpenFE execution tooling,
i.e. through the `quickrun method <https://docs.openfree.energy/en/latest/guide/execution/quickrun_execution.html>`_.
You can find additional information and examples on how to run simulations of the entire network in the "Running the simulations" section of our `CLI tutorial <https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html>`_.
In this study, we will use a slightly modified approach to our CLI Tutorial, allowing for execution of a single task per HPC job.

Here is an example of a very simple script that will create and submit a separate job script (`*.job` named file) for every alchemical transformation (for the simplest SLURM use case):

.. code-block:: bash

   for file in network_setup/transformations/*.json; do
     relpath="${file:30}"  # strip off "network_setup/"
     dirpath=${relpath%.*}  # strip off final ".json"
     jobpath="network_setup/transformations/${dirpath}.job"
     if [ -f "${jobpath}" ]; then
       echo "${jobpath} already exists"
       exit 1
     fi
     for repeat in {0..2}; do
       cmd="openfe quickrun ${file} -o results_${repeat}/${relpath} -d results_${repeat}/${dirpath}"
       echo -e "#!/usr/bin/env bash\n${cmd}" > "${jobpath}"
       sbatch "${jobpath}"
     done 
   done

Please reach out to the openfe team if you have any questions on how to adapt this script to your internal needs, we would be happy to assist with this.


.. _failed_edges:

Handling Failed Edges
=====================

It is possible that some of the simulations in the network fail.
More specifically, not all repeats of an edge may finish successfully, or an edge can fail entirely.
You should follow this strategy for dealing with those failures:

**1. Non-reproducible failures**

* When at least one repeat of the edge completed successfully
* Keep a log of the failure
* Rerun the failed job(s) up to 3 times.
* If the simulation repeat is still failing after 3 times, and the failure is due to NaN errors, then treat this as a **reproducible failure**.

**2. Reproducible failures**

* When all repeats of an edge failed
* Keep a log of the failure
* Do not rerun this edge
* If this is a redundant edge:

  * Remove this edge from the network and carry out the analysis without this edge

* If this is a non-redundant edge (meaning that removing this edge would lead to a disconnected graph):

  * Add a new edge to the network. We are currently working on a script that will automatically find a suitable new edge.
 
Identifying Failed Edges
------------------------

Failed edges can be identified in a few different ways:

1. Assuming you have no jobs running, there should be a 1-1 mapping between result ``json`` files and job work directories, for example:

.. code-block:: bash

   $ ls
   easy_rbfe_lig_ejm_31_complex_lig_ejm_46_complex/
   easy_rbfe_lig_ejm_31_complex_lig_ejm_46_complex.json
   easy_rbfe_lig_ejm_31_complex_lig_ejm_47_complex/
   easy_rbfe_lig_ejm_31_complex_lig_ejm_47_complex.json
   easy_rbfe_lig_ejm_31_complex_lig_ejm_48_complex/
   easy_rbfe_lig_ejm_31_complex_lig_ejm_48_complex.json
   easy_rbfe_lig_ejm_31_complex_lig_ejm_50_complex/
   easy_rbfe_lig_ejm_31_complex_lig_ejm_50_complex.json
   easy_rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent/
   easy_rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent.json

If a directory like ``easy_rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent/`` exists but the corresponding ``json`` file doesn't (``easy_rbfe_lig_ejm_31_solvent_lig_ejm_42_solvent.json``) then we know the edge failed.
The directory should be removed and the job should be resubmitted (depending on the failure type as discussed above).

2. The extract results script mentioned :ref:`here <inspecting results>` and the cleanup script mentioned :ref:`here <post-simulation cleanup>` both include output of which folder(s) and ``json`` file(s) contain errors and can be removed prior to starting new jobs.

Fixing broken networks
----------------------

If removing of a reproducibly failing edge leads to a disconnected graph, new edges need to be added in order to fix the broken network.
The fixing of disconnected networks is carried out using the script provided under
`utils/fix_networks.py <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/utils/fix_networks.py>`_.

Note that the script will throw an error if at least one repeat completed successfully. In that case (non-reproducible failure) we recommend re-running the failed jobs (see :ref:`Handling failed edges <failed_edges>`

Here an example of how to run the script:

.. code-block:: bash

   python fix_networks.py --input_alchem_network_file alchemicalNetwork/alchemical_network.json --result_files results_*/*json --output_extra_transformations new_edges

The script takes as inputs the `AlchemicalNetwork` from the original setup, the result `.json` files, and the folder name for storing the outputs.
This command will create a folder (named `new_edges` as specified using the `--output_extra_transformations` flag) that contains a `transformations` folder with a separate .json file for the solvent and complex legs for every new edge connecting the previously broken network.
The `new_edges` folder also contains a `ligand_network.graphml` file that is a serialized version of the full LigandNetwork (combining the old and the new `LigandNetwork` s) as well as an `alchemical_network.json` file containing the serialized version of the new `AlchemicalNetwork`.

Once the inputs for the edges to fix the broken network have been created, you can submit those calculations as described in the :ref:`Simulation execution section <simulation_execution>`.
Note that you will have to update the filepath to point to the new input .json files, e.g.

.. code-block:: bash

   for file in new_edges/transformations/*.json; do
     relpath="${file:30}"  # strip off "network_setup/"
     dirpath=${relpath%.*}  # strip off final ".json"
     jobpath="new_edges/transformations/${dirpath}.job"
     if [ -f "${jobpath}" ]; then
       echo "${jobpath} already exists"
       exit 1
     fi
     for repeat in {0..2}; do
       cmd="openfe quickrun ${file} -o results_${repeat}/${relpath} -d results_${repeat}/${dirpath}"
       echo -e "#!/usr/bin/env bash\n${cmd}" > "${jobpath}"
       sbatch "${jobpath}"
     done
   done


Inspecting Results
==================

.. _inspecting results:


.. note::
   A separate script will be provided for gathering relevant FE output data in Phase 3.


Due to the slightly modified simulation execution layout, using `openfe gather` will not work in openfe v1.0.1.

Instead we recommend using the scripts we bundle in this repository.

**Extracting ddG results**

You can extract the ddG results in the following manner from the directory containing `results_0`, `results_1`, and `results_2`.


.. code-block:: bash

   wget https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/extras/extract_results.py
   python extract_results.py


This will create a file named `ddg.tsv` from the individual simulation repeats.

If the simulations are complete, you can plot the `ddg.tsv` results
against the Schrodinger results downloaded from the relevant systems' ligand prediction file
from: https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/21_4_results/ligand_predictions

For example, for the JACS TYK2 set:


.. code-block:: bash

   wget https://raw.githubusercontent.com/schrodinger/public_binding_free_energy_benchmark/v2.0/21_4_results/ligand_predictions/jacs_set/tyk2_out.csv
   wget https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/extras/plot_results.py
   python plot_results.py --calculated ddg.tsv --experiment tyk2_out.csv


This will automatically generate plots of OpenFE vs experiment and FEP+ vs experiment using cinnabar.


Simulation Cleanup
==================

.. _post-simulation cleanup:

The post-simulation cleanup script will reduce the amount of data you need to store after your simulations.
It does not delete any data required for analysis.
To obtain the script, download it with:

.. code-block:: bash

   $ curl -LOJ https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/results_cleanup.py

The cleanup script requires the path to the ``result.json`` file and also accepts a list of ``result.json`` files.
Additionally, the script will detect if the result data has already been reduced and will skip processing that file.
This means that you can run the script repeatedly with a glob as results are finished.
For example:

.. code-block:: bash

   $ micromamba activate openfe-benchmark  # don't forget to activate your conda environment first!
   $ python results_cleanup.py path/to/results/*.json

.. warning::
   Do not move data around before running the script.
   The script assumes that the result directory has not been moved.

.. warning::
   Do not run multiple instances of the script at once using a wild card glob, i.e. ``*.json``.
   The script does not use file locks to ensure multiple instances are not operating on the same file.

To save space as simulations complete, consider adding a line at the end of your job submission script that runs the cleaning script.

.. code-block:: bash
   
   openfe quickrun "${file}" -o results_"${repeat}"/"${relpath}" -d results_"${repeat}"/"${dirpath}"
   python results_cleanup.py results_"${repeat}"/"${relpath}"

The cleanup script will delete:

* Redundant structural analysis data
* Redundant input data
* Simulation ``.nc`` & checkpoint ``.chk`` files 

.. note:: 
      
         The simulation ``.nc`` file is subsampled and for each lambda window, an ``.xtc`` trajectory file is created
   

Compute Requirements
====================

To run the benchmark simulations following **GPU hardware** will be required:

* An NVIDIA GPU (CUDA 10.2 or above compatible)
    * In non-exclusive compute mode
* A single GPU ID assigned per `openfe quickrun` execution
    * e.g. by setting `CUDA_VISIBLE_DEVICE` if necessary
* Estimated **standard** transformation compute time:
    * Approximately 8-12 GPU hours per complex transformation repeat
    * Approximately 1-2 GPU hours per solvent transformation repeat
* Estimated **net charge** transformation compute time:
    * Approximately 4-7 GPU days per complex transformation repeat
    * Approximately 8-12 GPU hours per solvent transformation repeat


Data Storage Requirements
=========================

.. _data storage requirements:

.. warning::
   You will need to keep any temporary data until the `post-simulation cleanup`_ script is made available.


**Temporary data:**

*Retention time:* until the `post-simulation cleanup`_ script is applied to your data.

*Estimated storage costs:*
  * Standard simulations: **5-10 GB** per triplicate repeat cycle
  * Net charge transformation: **40-80 GB** per triplicate repeat cycle

This data contains the full simulation outputs, including large netcdf trajectories with 1 ps snapshots of the coordinates and energies. A cleanup script will be provided to extract relevant long-term data from these outputs.

**Long term data:**

*Retention time:* until the benchmark manuscript is accepted for publication.
*Estimated storage costs:* **sub 500 mb** per triplicate repeat cycle.

This long term data will be extracted from above-mentioned temporary data using the `post-simulation cleanup`_ script.

It will include:
  * Reduced potential arrays for free energy analysis
  * OpenFE output JSON files
  * PNGs from structural analysis
  * XTC coordinate trajectories for each lambda window containing 21 evenly spaced frames

The data will be kept around by each partner for further processing should it be necessary as part of the manuscript writing process.


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

