import pandas as pd
from cinnabar import stats
import numpy as np
from typing import Optional

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
    statistics_out = []
    for statistic in statistics:
        s = stats.bootstrap_statistic(x,
                                      y,
                                      xerr,
                                      yerr,
                                      statistic=statistic,
                                      include_true_uncertainty=bootstrap_x_uncertainty,
                                      include_pred_uncertainty=bootstrap_y_uncertainty)
        statistics_out.append([statistic, s[statistic_type], s['low'], s['high']])
    return statistics_out

# Only calculating statistics for approprietly large systems
# According to Hahn et al. systems with a minimum dynamic range of 3 kcal/mol 
# and a minimum of 16 ligands (https://doi.org/10.33011/livecoms.4.1.1497)
min_dynamic_range = 3
min_number_ligands = 16

# Statistics to calculate for the datasets
statistics = ["RAE", "RMSE", "MUE", "R2", "rho", "KTAU"]

url_all = 'https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/processed_results/combined_pymbar3_calculated_dg_data.csv'
# Results after rerunning two of the Merck set/PFKFB3 edges after fixing a bug in the Kartograf atom mapper
url_rerun = 'https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/processed_results/reruns/rerun_pymbar3_calculated_dg_data.csv'

df = pd.read_csv(url_all)
df_rerun = pd.read_csv(url_rerun)
df_no_pfkfb3 = df[df["system name"] != 'pfkfb3']
df = pd.concat([df_no_pfkfb3, df_rerun])
df.reset_index(inplace=True, drop=True)

stats_openfe = {}
for s in df['system group'].unique():
    df_s = df[df['system group'] == s]
    df_s.reset_index(inplace=True, drop=True)
    stats_dict = {}
    for t in df_s['system name'].unique():
        df_t = df_s[df_s['system name'] == t]
        df_t.reset_index(inplace=True, drop=True)

        # Calculate the dynamic range of the dataset
        dynamic_range = abs(max(df_t['Exp DG (kcal/mol)']) - min(df_t['Exp DG (kcal/mol)']))
    
        # Do not calculate statics for systems with few ligands (specified in min_number_ligands)
        # or a small dynamic range (specified in min_dynamic_range)
        if dynamic_range < min_dynamic_range or len(df_t) < min_number_ligands:
            continue

        # Store statistics including the lower and upper 95% confidence interval        
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

print(stats_openfe)
