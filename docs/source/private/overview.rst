.. _private_overview:

**************************
Private Dataset Benchmarks
**************************


This page outlines the plans and instructions for the private dataset benchmarking portion of the 2024 OpenFE Industry Benchmark study.


Overview
********

The Private Dataset Benchmark focuses on validating the `OpenFE toolkit <https://docs.openfree.energy/en/stable/>`_
using blinded internal datasets from each industry partner. The aim of this study is to evaluate how well the OpenFE toolkit
is likely to perform in production cases.

Here each industry partner will be expected to:

* Gather suitable protein-ligand binding datasets with associated experimental data.
* Use the OpenFE tooling to estimate binding free energies on one or more internal dataset.
* Report blinded correlation plots of free energy estimates against experimental data.

Blinded data will be collated alongside any relevant meta analyses and published in an appropriate location.


Before you start
****************


Please ensure that you use the same environment used to run the public dataset OpenFE
benchmarks. Please see our benchmarks specific :ref:`benchmark openfe install instructions <installation>`
for more details.


Phase 1: Simulating Private Benchmark Sets
******************************************

Dataset Selection
=================

Each participating industry partner will be expected to select and prepare
benchmarks sets from their own internal projects.

Our aim is to try to assess how well OpenFE performs on the types of projects
that are being worked on in industry. You are therefore asked to pick **a minimum of one**
dataset for this task.

The datasets should:

* Have assay data for each ligand in the set.
* Represent the types of campaigns you run internally; this can either be fully internal or come from pre-published work.
* Not be part of any commonly published free energy benchmark sets.
* Where possible, follow `best practices in benchmark set selection <https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1497>`_.

You are also **encouraged to avoid** datasets which are known to not work with OpenFE.
This includes:

* Membrane-containing systems.
* Ligand series that undergo cyclisation.
* Ligand series that undergo scaffold hopping (i.e. very small to no conserved core).

Should you end up using a system with known challenges, we would ask you to describe the challenging nature of the dataset upon submission.


Dataset Preparation
===================

How the datasets are prepared is fully up to you. We encourage that you use lessons learnt from the
public set :ref:`input preparation step <input-preparation>` as a guide on how to prepare systems
for use with OpenFE. We expect that you will have a PDB for your protein, an SDF for your ligands,
and optionally another SDF for your cofactors. We encourage you to use the
:ref:`input validation script <input-validation>` to check that your inputs are ready for use with
OpenFE.


The only additional requirements for dataset preparations are:

* Any ligand names should be anonymised. The OpenFE data gathering scripts will not modify existing names, instead assuming them to be anonymous.
* If possible, you should record any methods used in dataset preparation, as suitable for the SI of a journal publication.


Running OpenFE Simulations
==========================

The same instructions as those from the :ref:`public datasets <public_phase2>` should be used here.

Lomap networks should be created using the script provided under
`utils/plan_rbfe_network.py <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/utils/plan_rbfe_network.py>`_.


.. code-block:: bash

   # If you don’t have cofactors
   python plan_rbfe_network.py --pdb protein.pdb --ligands ligands.sdf --output network_setup

   # If you have cofactors
   python plan_rbfe_network.py --pdb protein.pdb --ligands ligands.sdf --cofactors cofactors.sdf --output network_setup


You should then execute each transformation using the `quickrun method <https://docs.openfree.energy/en/latest/guide/execution/quickrun_execution.html>`_.
Below is an example script that will create and submit each job to a SLURM cluster scheduler:


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


My Simulations Are Failing, What Do I Do?
=========================================


Unfortunately for these benchmarks, the OpenFE team will only be able to provide help in a
limited manner (i.e. we will not be able to look at your structures).

To help you, we have created a :ref:`preparing and debugging simulations FAQ <prep_and_debug>` with
some common issues you may encounter.


Cleaning Results
================


.. note::
   Please keep all post-cleanup data around for analysis until the end of the benchmarks (i.e.
   after publication).


The OpenFE tools are known to generate a lot of data by default (something we are looking to fix!).

We recommend that folks use the :ref:`simulation cleanup <post-simulation cleanup>` script to clean up
unnecessary data.


Inspecting Results
==================


.. note::
   This section does not describe how data will be gathered by the OpenFE team for further analysis.
   A separate script will be provided for this purpose. See the :ref:`data gathering information <private-data-gather>`
   for more details.


If you wish to look at your results, you can use the `extract_results.py` script used in the
public dataset benchmarks:


.. code-block:: bash

   wget https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/extras/extract_results.py
   python extract_results.py


This will provide both dG and ddG outputs for you to further manipulate.

As we cannot tell what format your experimental results are in, we do not provide a plotting script at this time
and encourage you to use your own internal plotting tools.


You are encouraged to share early results with everyone on the #industry-benchmarking slack channel!


Handling Failed Edges
=====================


.. note::
   Please keep a note of any failed edges and report them when you submit results.


You should handle failed edges in the same way as the :ref:`public datasets <failed_edges>`.

Reproducible failures which result in broken networks can be fixed using the same method as the :ref:`public datasets <failed_edges>`.


Phase 2: Data Gathering
***********************

What data will we gather?
=========================

Like during the public phase of this benchmarking project we need to gather some general data about the alchemical network you have run. More information
about the script which extracts the data and the contents of the final data package to be shared with the OpenFE team can
be found :ref:`here <gathering_of_results>`.

As this phase involves private datasets we also require you to disclose the experimental data used to validate the accuracy
of the simulations in a specific CSV format which is shown here for some example data:

.. csv-table:: Example Data
   :file: experimental_data.csv
   :header-rows: 1

.. warning::
    The names of the ligands must match those used when planning your alchemical network and
    therefore those in the gathered results. If the names are private and you hid them using the
    :ref:`data gathering script <gathering_of_results>`, an additional :ref:`script <rename-csv-ligands>` is provided which can translate
    those in the experimental CSV to match the anonymised names.

Notes on generating the CSV data:

* Experimental assay values should be supplied as ``nanomolar affinities``.
* If experimental error is not available a value of ``-1`` should be used.
* The annotation column should be used to note anything different about this ligand, such as the assay changing or the measurement being near the limit.


.. _rename-csv-ligands:

Renaming ligands
~~~~~~~~~~~~~~~~
The ``rename_exp_data.py`` script can be used to hide any private ligand names in your created experimental data file
using the name mapping which was created by the ``data_gathering.py`` script.

.. code-block:: bash

   wget https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/rename_exp_data.py
   python rename_exp_data.py --experimental-data my_exp_data.csv -name-mapping-file ligand_name_mapping_PRIVATE.json --output blinded_exp_data.csv

Private datasets info form
~~~~~~~~~~~~~~~~~~~~~~~~~~
We also wish to gather additional information about the private datasets:

* Details on system preparation
* Estimated benchmark difficulty (i.e. "easy", "medium", "hard")
* Additional details which may impact simulation difficulty, e.g. "likely water sampling issues" or "ions in the binding site".
* Broad assay description (e.g. "Kd from ITC") and additional details on the experimental assay.
* Compute hardware description, e.g. what type of GPU was used.
* Any information you have on systematic edge failures

We provide a template text file for you to fill in all the information. You can download this file using:

.. code-block:: bash

   wget https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/main/industry_benchmarks/utils/info_form_private_sets.txt

You can open and edit the file using e.g. VIM or Microsoft Word. Please fill out this form for each private dataset and upload
the file (as a pdf or text file) together with the other results in your zenodo upload.


How will this be gathered?
==========================

.. _private-data-gather:

We anticipating gathering the data through:

1. A script provided by the OpenFE team that will extract relevant information from your simulations.
2. A text file for you to fill in with information.
3. A CSV file with experimental data that you will need to prepare.


The text form and data gathering scripts are currently being prepared by the OpenFE team. Please come back soon for more details!
