## System preparation details

* System capped with pymol 2.5.4
* Move the TER card after the NME
* Rename the GLY HA atoms on the first residue to HA3 and HA2 respectively
* Manually rename the HD2 atom on the terminal HIS
* Removed ligand `30_protonated` per:
```
https://static-content.springer.com/esm/art%3A10.1038%2Fs42004-023-01019-9/MediaObjects/42004_2023_1019_MOESM2_ESM.pdf
"Ligand 30 has a pyridine ring that was modeled as being neutral. EpiK predicts the pKa of the piperidine nitrogen to be 6.45, implying that, depending on the pH of the assay, the protonated form could be a significant population in solvent. The pH of the assay was not reported, but BACE1 is primarily active between pHs of 4 and 555. When only the neutral form was modeled in the map, the binding free energy of ligand 30 was predicted to be too negative by 1.7 kcal mol−1 . Adding the protonated form of ligand 30 and applying the pKa"
```
