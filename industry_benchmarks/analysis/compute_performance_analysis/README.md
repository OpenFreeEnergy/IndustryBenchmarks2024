## Compute performance analysis

This folder contains CSV data for the ns/day performance of the public dataset using representative edges.
Code and inputs for running this analysis can be found here: https://github.com/OpenFreeEnergy/performance_benchmarks

## Analysis design

The edge with the largest amount of atoms being transformed in each dataset was chosen as a representative edge.
A short 500 ps simulation was run to obtain a stable ns/day value for each edge.

All edges were run on the same hardware configuration using the OpenFE v1.0.1 conda-lock file, with OpenMM 8.1.1 and cuda 11.8.

## Datasets

- `L40S_performance_analysis.csv`: Performance on a L40S GPU.
