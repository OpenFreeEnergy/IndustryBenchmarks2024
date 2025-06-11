import pandas as pd
from cinnabar import stats
import numpy as np
from typing import Optional


def _filter_ross_df(df_ross, df_openfe):
    # Filter out systems and ligand variations that were not run by OpenFE
    indxs = []
    for inx, row in df_ross.iterrows():
        df_s = df_openfe[df_openfe['system group'] == row['system group']]
        if row['system name'] not in df_s['system name'].unique():
            indxs.append(inx)
        if row['ligand name'] not in df_s['ligand name'].unique():
            indxs.append(inx)
    df_ross_filtered = df_ross.drop(indxs)
    df_ross_filtered.reset_index(inplace=True, drop=True)
    return df_ross_filtered


def _get_statistics(
    x: np.ndarray,
    y: np.ndarray,
    xerr: Optional[np.ndarray] = None,
    yerr: Optional[np.ndarray] = None,
    statistics: list = ["RMSE", "MUE"],
    bootstrap_x_uncertainty: bool = False,
    bootstrap_y_uncertainty: bool = False,
    statistic_type: str = "mle",
):
    statistics_out = {}
    for statistic in statistics:
        s = stats.bootstrap_statistic(x,
                                      y,
                                      xerr,
                                      yerr,
                                      statistic=statistic,
                                      include_true_uncertainty=bootstrap_x_uncertainty,
                                      include_pred_uncertainty=bootstrap_y_uncertainty)
        statistics_out[statistic] = s[statistic_type]
        statistics_out[f'{statistic}_low'] = abs(s[statistic_type] - s['low'])
        statistics_out[f'{statistic}_high'] = abs(s[statistic_type] - s['high'])
    return statistics_out


def _get_subset_dataframe(df, column_name, value):
    df_sub = df[df[column_name] == value]
    df_sub.reset_index(inplace=True, drop=True)
    return df_sub


def _get_statistics_dict(
    df: pd.DataFrame,
    statistics: list[str] = ["RAE", "RMSE", "MUE", "R2", "rho", "KTAU"],
    min_dynamic_range: int = 3,
    min_number_ligands: int = 16,
):
    """
    Function to calculate statistics for all datasets.
    Only calculating statistics for appropriately large systems
    According to Hahn et al. systems with a minimum dynamic range of 3 kcal/mol
    and a minimum of 16 ligands (https://doi.org/10.33011/livecoms.4.1.1497)

    Parameters
    ----------
    df: pd.DataFrame
      DataFrame including data
    statistics: list[str]
      Statistics to calculate for the datasets
    min_dynamic_range: int
      minimum dynamic range (in kcal/mol)
    min_number_ligands: int
      minimum number of ligands per dataset

    Returns
    -------
    results: dict[str, dict[str, dict[str, float]]]
      Results dictionary of the different sets and individual datasets within each set.
      Contains the name of the static and a list of the value for that statistic,
      the lower, and the upper confidence interval.

    """

    stats_openfe = {}
    for s in df['system group'].unique():
        df_s = _get_subset_dataframe(df, 'system group', s)
        stats_dict = {}
        for t in df_s['system name'].unique():
            df_t = _get_subset_dataframe(df_s, 'system name', t)
    
            # Calculate the dynamic range of the dataset
            dynamic_range = abs(max(df_t['Exp DG (kcal/mol)']) - min(df_t['Exp DG (kcal/mol)']))
        
            # Do not calculate statics for systems with few ligands (specified in min_number_ligands)
            # or a small dynamic range (specified in min_dynamic_range)
            if dynamic_range < min_dynamic_range or len(df_t) < min_number_ligands:
                continue
    
            # Store statistics including the lower and higher 95% confidence interval        
            stats_out = _get_statistics(
                df_t['Exp DG (kcal/mol)'], 
                df_t['DG (kcal/mol)'], 
                df_t['Exp dDG (kcal/mol)'], 
                df_t['uncertainty (kcal/mol)'],
                statistics=statistics,
            )
            
            stats_dict[t] = stats_out 
        if stats_dict:
            stats_openfe[s] = stats_dict

    return stats_openfe


url_all = 'https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/processed_results/combined_pymbar3_calculated_dg_data.csv'
# Results after rerunning two of the Merck set/PFKFB3 edges after fixing a bug in the Kartograf atom mapper
url_rerun = 'https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/processed_results/reruns/rerun_pymbar3_calculated_dg_data.csv'
# FEP+ results from Ross et al.
url_ross = 'https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/schrodinger_21_4_results/combined_schrodinger_dg.csv'

# Get the DataFrame for the OpenFE results
df = pd.read_csv(url_all)
df_rerun = pd.read_csv(url_rerun)
df_no_pfkfb3 = df[df["system name"] != 'pfkfb3']
df = pd.concat([df_no_pfkfb3, df_rerun])
df.reset_index(inplace=True, drop=True)

# Get the DataFrame for the Ross et al (FEP+) results
df_ross = pd.read_csv(url_ross)
# Rename some columns to match OpenFE csv file
df_ross = df_ross.rename(columns={
    "Ligand name": "ligand name",
    "Exp. dG (kcal/mol)": "Exp DG (kcal/mol)",
    "Pred. dG (kcal/mol)": "DG (kcal/mol)",
    "Exp. dG error (kcal/mol)": "Exp dDG (kcal/mol)",
    "Pred. dG std. error (kcal/mol)": "uncertainty (kcal/mol)",
})
# Filter the Ross Dataframe to only include ligands run by OpenFE
df_ross_filtered = _filter_ross_df(df_ross, df)

# Get the statistics for all systems
statistic_openfe = _get_statistics_dict(df)
statistic_ross = _get_statistics_dict(df_ross_filtered)
