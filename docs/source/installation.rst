.. _installation:

OpenFE Installation
*******************

The page has information for installing ``openfe``, and testing that your ``openfe`` installation is working.

For this benchmarking study, we would like to ask that all participants use a specific installation of OpenFE (v1.0.1) alongside selected pinned tooling dependencies (e.g. openff-toolkit v0.15.2, openmm v8.1.1, openmmforcefields v0.13.0, openmmtools v0.23.1).

To do this, we ask participants to install a version of openfe with a pre-defined environment. This means that you may likely have to re-install any existing `openfe` installs. This can be achieved in two different ways, using a conda lock file or a single file installer as described below.

Please reach out to the OpenFE team should installing a new openfe version not be easily feasible or if you need help with your installation.

.. note::

   We recommend using a ``conda-lock`` file unless you are unable to download packages.



``conda-lock`` file
===================

.. _conda-lock: https://github.com/conda/conda-lock?tab=readme-ov-file#conda-lock

A `conda-lock`_ file is a cross platform way of specifying a conda environment that specifies packages in a reproducible way.
Unlike the single file installer, an internet connection is required to install from a ``conda-lock`` file.
We recommend the use of a ``conda-lock`` file when the same conda environment is required across different systems.

See https://github.com/conda/conda-lock?tab=readme-ov-file#conda-lock for more information on ``conda-lock``.

The `conda-lock` file for OpenFE version v1.0.1 can be downloaded with ::

  $ curl -LOJ https://github.com/OpenFreeEnergy/openfe/releases/download/v1.0.1/conda-lock-openfe-1.0.1.yml

We recommend using ``micromamba`` when working with ``conda-lock`` files ::

  $ micromamba create -n openfe --file openfe-1.0.1-conda-lock.yml
  $ micromamba activate openfe

If ``micromamba`` is not available, ``conda`` or ``mamba`` may be used after installing ``conda-lock``.

.. note::

   You will likely need to install ``conda-lock``.
   We recommend installing ``conda-lock`` in a new virtual environment using either `conda` or `mamba`.
   This will reduce the chance of dependency conflicts ::

       $ # Install conda lock into a virtual environment
       $ conda create -n conda-lock -c conda-forge conda-lock
       $ # Activate the environment to use the conda-lock command
       $ conda activate conda-lock

Create a conda environment from the lock file and activate it::

  $ conda-lock install -n openfe conda-lock-openfe-1.0.1.yml
  $ conda activate openfe

For additional details, please visit the `Installation page <https://docs.openfree.energy/en/latest/installation.html>`_ in the OpenFE documentation.

Single file installer
=====================

.. note::

   We recommend using a ``conda-lock`` file unless you are unable to download packages.
   The error may look something like this ::

       ERROR:root:Retrying (Retry(total=2, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NameResolutionError("<urllib3.connection.HTTPSConnection object at 0x7bc5c3e75670>: Failed to resolve 'conda.anaconda.org' ([Errno -2] Name or service not known)")'

   Only then should you use the single file installer.
   If ``mamba`` or ``micromamba`` are able to install packages from ``conda-forge`` then we recommend using the ``conda-lock`` file.


.. warning::

   The single file installer may modify your ``.bashrc`` in a way that requires manual intervention to access your previous ``conda`` installation 

.. _releases on GitHub: https://github.com/OpenFreeEnergy/openfe/releases

Single file installers are available for x86_64 Linux and MacOS.
They are attached to our `releases on GitHub`_ and can be downloaded with a browser or ``curl`` (or similar tool).
For example, the Linux installer can be downloaded with ::

  $ curl -LOJ https://github.com/OpenFreeEnergy/openfe/releases/download/v1.0.1/OpenFEforge-1.0.1-Linux-x86_64.sh

And the MacOS (x86_64) installer ::

  $ curl -LOJ https://github.com/OpenFreeEnergy/openfe/releases/download/v1.0.1/OpenFEforge-MacOSX-x86_64.sh

And the MacOS (arm64) installer ::

  $ curl -LOJ https://github.com/OpenFreeEnergy/openfe/releases/download/v1.0.1/OpenFEforge-MacOSX-arm64.sh

The single file installer contains all of the dependencies required for ``openfe`` and does not require internet access to use.

Both ``conda`` and ``mamba`` are also available in the environment created by the single file installer and can be used to install additional packages.
The installer can be installed in batch mode or interactively  ::

  $ chmod +x ./OpenFEforge-1.0.1-Linux-x86_64.sh # Make installer executable
  $ ./OpenFEforge-1.0.1-Linux-x86_64.sh # Run the installer

After the installer completes, close and reopen your shell.
To check if your path is setup correctly, run ``which python`` your output should look something like this ::

   (base) $ which python
   /home/mmh/openfeforge/bin/python

.. note::
   Your path will be different, but the important part is ``openfeforge/bin/python``

For additional details, please visit the `Installation page <https://docs.openfree.energy/en/latest/installation.html>`_ in the OpenFE documentation.


Testing your installation
=========================

To make sure everything is working, run the tests ::

  $ pytest --pyargs openfe openfecli

The test suite contains several hundred individual tests. This will take a
few minutes, and all tests should complete with status either passed,
skipped, or xfailed (expected fail).

With that, you should be ready to use ``openfe``!
