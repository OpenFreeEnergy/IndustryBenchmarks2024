.. _contributing-inputs:

Contributing inputs
###################

This page will walk you through the process of working with this Github repository and uploading your prepared input files.

**1. Clone the github repository**

.. code-block:: bash

   git clone https://github.com/OpenFreeEnergy/IndustryBenchmarks2024.git


**2. Create and work off of a remote branch**

.. code-block:: bash

   git checkout -b my_remote_branch

Where ``my_remote_branch`` is the name for you remote branch, e.g. ``prepare_systems_tyk2_hif2a``.

Now you can add your prepared files under the ``industry_benchmarks/inputs/prepared_structures`` subfolder.

**3. Push the prepared input files to the Github repository**

Once the preparation is completed, you can upload the files onto Github:

.. code-block:: bash

   # Add all your files
   git add <your file 1>  <your file 2> <your file ...>
   # Create a commit with a meaningful commit message
   git commit -m 'Prepared input files for system X'
   # Push the commit to GitHub
   git push --set-upstream origin <my_remote_branch>


**4. Create a Pull Request**

In a next step, you can create a Pull Request on Github.

ToDo:
      Add pictures of how it would look like

      Add link to PR template for input file submission



