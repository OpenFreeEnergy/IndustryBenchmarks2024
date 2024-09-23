.. _prep_and_debug:

***********************************************
FAQ: Preparing and Debugging OpenFE Simulations
***********************************************

This FAQ contains a list of issues you may encounter when preparing and
running simulations with OpenFE. This list will be updated as feedback
is provided by benchmark partners.


1. Ligand Alignment
*******************

Ligands in a series should be pre-aligned when pass to OpenFE tools.
This means that regions of the ligand that should be part of a common
core, should have the same coordinates.

In these benchmarks we will use the Kartograf atom mapper. Any ligand
regions which do not align in 3D space will be unmapped and treated
as unique atoms in the alchemical topology.


2. Protein Structures Ready for Simulation
******************************************

OpenFE does not run a very long equilibration during the RBFE calculation.
As such, a protein structure with clashes will likely fail with the Protocol
reporting "NaN" errors.

You should make sure that your protein structure is ready for MD simulation
before passing it to OpenFE.


3. All My Complex Edges Fail But Not The Ligands
************************************************

This indicates a likely issue with the protein structure, either a clash
with the ligand, or a clash within the protein structure itself.

To check this we suggest running a conventional MD simulation and seeing if the system remains stable.
You can use the `OpenFE MD Protocol for this purpose <https://docs.openfree.energy/en/stable/tutorials/md_tutorial.html>`.


4. My Ligand Transformations Are Failing
****************************************

This likely indicates a bad mapping.

We would encourage you to `visualize your mappings <https://docs.openfree.energy/en/stable/cookbook/ligandnetwork_vis.html>`.
If you see unusual atom assignments, please get in touch with the OpenFE team and we will guide you through how to assign alternative mappings.
