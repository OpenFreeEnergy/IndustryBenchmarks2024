.. _contributing-inputs:

*******************
Contributing inputs
*******************

This page will walk you through the process of contributing your prepared input files.

Checklist
*********

When submitting input files, you should have completed the following:

* created a directory for each system based on the `template directory <https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/input_structures/prepared_structures/template>`_
* named and placed the directory with the following pattern: `input_structures/prepared_structures/<set_name>/<system_name>`

For each submitted directory:

* filled in system preparation details in the `PREPARATION_DETAILS.md` file
* added a PDB file with the protein named `protein.pdb`
* removed any cofactors from the PDB file and placed them in `cofactors.sdf`
* kept cystallographic waters and metals in the PDB file
* tested the PDB and, if available, cofactors.sdf file using the [input validation script][input_validation]
* copied over the ligands SDF file to a file named `ligands.sdf`
* copied over the edges CSV file to a file named `edges.csv`


Detailed Instructions
*********************

**1. Fork the GitHub repository**

First, you can create a fork of the GitHub repository. You will then make the changes (adding the files) first in the forked repository and afterwards propose the changes to the upstream repository (`OpenFreeEnergy/IndustryBenchmarks2024`).

* In the `OpenFreeEnergy/IndustryBenchmarks2024` repository in the top-right corner of the page, click **Fork**.
* Under "Owner", select an owner for the forked repository (likely yourself). The repository name can stay as given, `IndustryBenchmarks2024`.
* Click **Create Fork**.

**2. Clone the forked GitHub repository**

You can clone the fork created in step 1, meaning making a local copy of the repository, to make it easier to add new files and push changes back to the remote repository on GitHub.com.

.. code-block:: bash

   git clone https://github.com/Owner/IndustryBenchmarks2024.git

where ``Owner`` is the name of the owner that was selected in step 1.

**3. Create and work off of a remote branch**

You can create a remote branch that will allow you to make changes to the repository (e.g. adding new files) without the changes immediately affecting the main reporitory.

.. code-block:: bash

   git checkout -b my_remote_branch

Where ``my_remote_branch`` is the name for you remote branch, e.g. ``prepare_systems_tyk2_hif2a``.

Now you can add your prepared files under the ``industry_benchmarks/input_structures/prepared_structures`` subfolder.

**4. Push the prepared input files to the Github repository**

Once the preparation is completed, you can upload the files onto GitHub.com:

.. code-block:: bash

   # Add all your files
   git add <your file 1>  <your file 2> <your file ...>
   # Create a commit with a meaningful commit message
   git commit -m 'Prepared input files for system X'
   # Push the commit to GitHub
   git push --set-upstream origin <my_remote_branch>

All files will be deposited within the directory ``input_structures/prepared_structures`` following the directory structure:
``input_structures/prepared_structures/<set_name>/<system_name>``.

Following files will need to be added to the respective directory for each system (based on the ``template directory``):

* ``PREPARATION_DETAILS.md``: Filled in system preparation details
* ``protein.pdb``: PDB file with the protein, including crystallographic waters and metals
* ``cofactors.sdf``: Cofactors that were moved from the original PDB file to this file
* ``ligands.sdf``: Ligand SDF
* ``edges.csv``: Edges CSV file

**5. Create a Pull Request into the upstream repository**

In a next step, you can create a Pull Request on GitHub. A Pull Request (PR) is a proposal to merge the changes from your remote branch into another branch, e.g. into the main codebase.
We created a PR template for you where you can tick off the checklist to ensure that all steps have been completed.
To create the PR:

* Go to your fork of the GitHub repository. After pushing the remote branch (step 3) there should now be a note that a new branch (with the name you had given it) has recently been pushed.
* Click on the green button that says "Compare & pull request". 
* In the "base repository" dropdown menu, select the upstream repository (``OpenFreeEnergy/IndustryBenchmarks2024``) and the "base branch" ``main``.
* Create the PR from the PR template provided, giving the PR a meaningful title and description.

