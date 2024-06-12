<!--
Thank you for opening a contribution to the OpenFE 2024 industry benchmark repository.
Below are a few things we ask you to kindly fill in and self-check before we
can accept your contribution. Please ignore any irrelevant sections.
-->

## Input submission checklist

<!--
If you are submitting prepared input files please indicate if you have
done the following:
-->

* [ ] created a directory for each system based on the [template directory][templates]
* [ ] named and placed the directory with the following pattern: `input_structures/prepared_structures/<set_name>/<system_name>`

For each submitted directory:

* [ ] filled in system preparation details in the `PREPARATION_DETAILS.md` file
* [ ] added a PDB file with the protein named `protein.pdb`
* [ ] removed any cofactors from the PDB file and placed them in `cofactors.sdf`
* [ ] kept cystallographic waters and metals in the PDB file
* [ ] tested the PDB and, if available, cofactors.sdf file using the [input validation script][input_validation]
* [ ] copied over the ligands SDF file to a file named `ligands.sdf` and removed any duplicate binding modes or protonation states, keeping the more favorable state

<!--
Here please add a summary of what changes you have made
-->
## Summary of changes

<!--
Also please indicate that you are happy to release these materials under the
combined MIT and CCBY-SA 4.0 licenses of this repository
-->
## Licensing agreement
* [ ] I agree to release these inputs under the combined MIT and CCBY SA 4.0 licenses of this repository

[templates]: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks
[input_validation]: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/tree/main/industry_benchmarks/utils/input_validation.py
