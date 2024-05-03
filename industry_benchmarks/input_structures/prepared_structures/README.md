# Remediated Input Structures

This directory holds remediated benchmark inputs ready for use with the OpenFE toolkit.

## Structure

Data is organized in the following structure:

```bash
  -- <benchmark_set>
      -- <benchmark_system>
          -- protein.pdb # PDB of protein + cocrystalized waters & ions
          -- ligands.sdf # SDF containing ligands to be transformed
          -- cofactors.sdf # (Optional) SDF containing any system cofactors
          -- edges.csv # CSV file containing the initial FEP+ study perturbation network
```

## Template

A [template directory](./template) is provided demonstrating how we anticipate what each sub-directory may contain.

