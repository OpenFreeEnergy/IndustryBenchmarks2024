# Dataset Analysis

The scripts provided in this folder should be used to process the industry benchmark datasets uploaded to Zenodo.

## Creating the analysis environments

You should use the recommended openfe-1.0.1 environment to run the analysis, this can be installed following the 
[guide](https://industrybenchmarks2024--157.org.readthedocs.build/en/157/installation.html#openfe-installation).


## Running the analysis

The scripts should be executed in the following order.

``1_download_and_extract_data.py``: This script will download the results archive from Zenodo and extract edge scores 
and DG estimates for each repeat and phase and collect them into a CSV file (`pymbar3_edge_data.csv`) along with the experimental data. It will
also calculate the cumulative DG estimates for each edge and report these into a separate CSV file 
(`pymbar3_cumulative_data.csv`). 

The script will automatically find the experimental data for private sets, for public sets you can provide the data as:

```bash
python 1_download_and_extract_data.py -z https://zenodo.org/records/14229113   \ 
                                      -p Merck -o merck_hif2a                  \
                                      -e  hif2a_automap_symbmcorr_out.csv
```