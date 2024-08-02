from pathlib import Path
import csv
import click
from cinnabar import Measurement, ReferenceState, FEMap
from cinnabar import plotting as cinnabar_plotting
from openff.units import unit
import numpy as np


def get_exp_data(filename: Path) -> dict[str, dict[str, float]]:
    """
    Fetch the experimental data from a Schrodinger Benchmark ligand file.
    
    Parameters
    ----------
    filename : pathlib.Path
      Pathlib object for the file with the experimental data.
    
    Returns
    -------
    experimental_data : dict[str, dict[str, float]]
      Dictionary of experimental data dictionaries for each ligand.
    """
    experimental_data = {}

    with open(filename, 'r') as fd:
        rd = csv.reader(fd, delimiter=',', quotechar='"')
        
        # Get the data indices
        headers = next(rd)
        name_idx = headers.index('Ligand name')
        exp_dg_idx = headers.index('Exp. dG (kcal/mol)')
        try:
            exp_dg_error_idx = headers.index('Exp. dG error (kcal/mol)')
        except ValueError:
            exp_dg_error_idx = None
        fep_dg_idx = headers.index('Pred. dG (kcal/mol)')
        fep_dg_error_idx = headers.index('Pred. dG std. error (kcal/mol)')
        
        # Loop through each ligand entry
        for row in rd:
            experimental_data[row[name_idx]] = {}
            experimental_data[row[name_idx]]['exp_dG'] = float(row[exp_dg_idx])
            if exp_dg_error_idx is not None:
                experimental_data[row[name_idx]]['exp_dG_err'] = float(row[exp_dg_error_idx])
            else:
                experimental_data[row[name_idx]]['exp_dG_err'] = 0
            experimental_data[row[name_idx]]['fep_dG'] = float(row[fep_dg_idx])
            experimental_data[row[name_idx]]['fep_dG_err'] = float(row[fep_dg_error_idx])

    return experimental_data
    

def get_calc_data(filename: Path) -> dict[str, dict[str, float]]:
    """
    Fetch the calculated data from an OpenFE ddG TSV file.
    
    Parameters
    ----------
    filename : pathlib.Path
      Pathlib object for the file with the calculated data.
    
    Returns
    -------
    calculated_data : dict[str, dict[str, float]]
      Dictionary of calculated data dictionaries for each ligand.
    """
    calculated_data = {}


    with open(filename, 'r') as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        headers = next(rd)
        for row in rd:
            tag = row[0] + "->" + row[1]
            calculated_data[tag] = {}
            calculated_data[tag]['ligand_i'] = row[0]
            calculated_data[tag]['ligand_j'] = row[1]
            calculated_data[tag]['ddG'] = float(row[2])
            calculated_data[tag]['ddG_err'] = float(row[3])

    return calculated_data


def get_femap(
    exp_data: dict[str, dict[str, float]],
    calc_data: dict[str, dict[str, float]]
) -> FEMap:
    """
    Create a cinnabar FEMap
    
    Parameters
    ----------
    exp_data: dict[str, dict[str, float]]
      The experimental data.
    calc_data: dict[str, dict[str, float]]
      The calculated data.
    
    Returns
    -------
    femap : cinnabar.FEMap
      A cinnabar FEMap object with generated absolute values.
    """
    # define a ground data
    ground = ReferenceState()

    fe_results = {'Experimental': {}, 'Calculated': []}

    # Lead the Measurements
    for entry in exp_data:
        m = Measurement(
            labelA=ground,
            labelB=entry,
            DG=exp_data[entry]['exp_dG'] * unit.kilocalorie_per_mole,
            uncertainty=exp_data[entry]['exp_dG_err'] * unit.kilocalorie_per_mole,
            computational=False
        )
        fe_results['Experimental'][m.labelB] = m

    for entry in calc_data:
        m = Measurement(
            labelA=calc_data[entry]['ligand_i'],
            labelB=calc_data[entry]['ligand_j'],
            DG=calc_data[entry]['ddG'] * unit.kilocalorie_per_mole,
            uncertainty=calc_data[entry]['ddG_err'] * unit.kilocalorie_per_mole,
            computational=True
        )
        fe_results['Calculated'].append(m)
    
    # Feed into the FEMap object
    femap = FEMap()

    for entry in fe_results['Experimental'].values():
        femap.add_measurement(entry)

    for entry in fe_results['Calculated']:
        femap.add_measurement(entry)

    femap.generate_absolute_values()

    return femap


def plot_femap(
    femap: FEMap,
    exp_data: dict[str, dict[str, float]],
    ddg_filename : str,
    dg_filename : str
 ) -> None:
    """
    Helper method to plot out ddG and dG plots.
    
    Parameters
    ----------
    femap : FEMap
      The cinnabar FEMap to plot
    exp_data : dict[str, dict[str, float]]
      Experimental data to use for shifting the absolute values
    ddg_filename : str
      The name of the ddg plot file.
    dg_filename : str
      The name of the dg plot file.
    """
    cinnabar_plotting.plot_DDGs(
        femap.to_legacy_graph(),
        figsize=5,
        filename=ddg_filename,
        xy_lim=[-5, 5],
    )
    
    shift = sum([i['exp_dG'] for i in exp_data.values()]) / len(exp_data)
    
    cinnabar_plotting.plot_DGs(
        femap.to_legacy_graph(),
        figsize=5,
        shift=shift,
        filename=dg_filename,
        xy_lim=[-15, -5],
    )
    

def plot_schrodinger_comparison(
    exp_data: dict[str, dict[str, float]],
    filename: str,
) -> None:
    """
    Helper method to plot out a dG comparison between experimental,
    and schrodinger.
    
    Parameters
    ----------
    exp_data: dict[str, dict[str, float]]
      The data loaded from the exp data file.
    filename : str
      The name of the plot file.
    """
    exp = []
    exp_err = []
    fep = []
    fep_err = []
    for entry in exp_data:
        exp.append(exp_data[entry]['exp_dG'])
        exp_err.append(exp_data[entry]['exp_dG_err'])
        fep.append(exp_data[entry]['fep_dG'])
        fep_err.append(exp_data[entry]['fep_dG_err'])
    
    shift = sum([i['exp_dG'] for i in exp_data.values()]) / len(exp_data)
    
    cinnabar_plotting._master_plot(
        np.asarray(exp),
        np.asarray(fep),
        xerr=np.asarray(exp_err),
        yerr=np.asarray(fep_err),
        xlabel='experimental',
        ylabel='fep+',
        statistics=["RMSE", "MUE", "R2", "rho"],
        title='Experiment vs FEP+',
        filename=filename,
        quantity= r"$\Delta$ G",
        bootstrap_x_uncertainty=False,
        bootstrap_y_uncertainty=False,
        statistic_type='mle',
        figsize=5,
        shift=shift,
        xy_lim=[-15, -5],
    )


@click.command
@click.option(
    '--calculated',
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    default=Path("ddg.tsv"),
    required=True,
    help=("Path to file output from extract_results.py"),
)
@click.option(
    '--experiment',
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    default=Path("experiment.csv"),
    required=True,
    help=("Path to Schrodinger Public Benchmarks ligands predictions file"),
)
@click.option(
    '--ddg_plot_filename',
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    default=Path("openfe_experiment_ddg.png"),
    required=True,
    help=("name of output openfe vs experiment ddG plot file"),
)
@click.option(
    '--dg_plot_filename',
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    default=Path("openfe_experiment_dg.png"),
    required=True,
    help=("name of output openfe vs experiment dG plot file"),
)
@click.option(
    '--dg_fepplus_plot_filename',
    type=click.Path(dir_okay=False, file_okay=True, path_type=Path),
    default=Path("schrodinger_experiment_dg.png"),
    required=True,
    help=("name of output schrodinger vs experiment dG plot file"),
)
def run(
    calculated : Path,
    experiment : Path,
    ddg_plot_filename : Path,
    dg_plot_filename : Path,
    dg_fepplus_plot_filename : Path
) -> None:
    """
    Get all the relevant plots.
    
    Parameters
    ----------
    calculated : Path
      Path to file output from extract_results.py.
    experiment : Path
      Path to Schrodinger Public Benchmarks ligands predictions file.
    ddg_plot_filename : Path
      Name of output openfe vs experiment ddG plot file.
    dg_plot_filename : Path
      Name of output openfe vs experiment dG plot file.
    dg_fepplus_plot_filename : Path
      Name of output schrodinger vs experiment dG plot file.
    """
    exp_data = get_exp_data(experiment)
    calc_data = get_calc_data(calculated)
    femap = get_femap(exp_data, calc_data)
    plot_femap(femap, exp_data, ddg_plot_filename, dg_plot_filename)
    plot_schrodinger_comparison(exp_data, dg_fepplus_plot_filename)


if __name__ == "__main__":
    run()
