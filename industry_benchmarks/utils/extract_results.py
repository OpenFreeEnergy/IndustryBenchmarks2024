import numpy as np
import glob
import json
import pathlib
import click


def get_names(result) -> tuple[str, str]:
    # Result to tuple of ligand names
    nm = list(result['unit_results'].values())[0]['name']
    toks = nm.split()
    if toks[2] == 'repeat':
        return toks[0], toks[1]
    else:
        return toks[0], toks[2]


@click.command
@click.option(
    '--results_0',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('results_0'),
    required=True,
    help=("Path to the directory that contains all result json files " 
          "for repeat 0, default: results_0."),
)
@click.option(
    '--results_1',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('results_1'),
    required=True,
    help=("Path to the directory that contains all result json files " 
          "for repeat 1, default: results_1."),
)
@click.option(
    '--results_2',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('results_2'),
    required=True,
    help=("Path to the directory that contains all result json files " 
          "for repeat 2, default: results_2."),
)
@click.option(
    '--output',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path('ddg.csv'),
    required=True,
    help="Path to the output csv file in which the DDG values will be stored "
         "default: ddg.csv.",
)
def extract(results_0, results_1, results_2, output):
    files_0 = glob.glob(f"{results_0}/*.json")
    files_1 = glob.glob(f"{results_1}/*.json")
    files_2 = glob.glob(f"{results_2}/*.json")
    if len(files_0) != len(files_1) != len(files_2):
        print(
            f'CAVE: Some edges did not finish: repeat 0: {len(files_0)} '
            f'edges, '
            f'repeat 1: {len(files_1)} edges, repeat 2: {len(files_2)} edges')
    edges_dict = dict()
    for file in files_0:
        with open(file) as stream:
            json_0 = json.load(stream)
        runtype = file.split('/')[1].split('_')[-2]
        molA, molB = get_names(json_0)
        edge_name = f'edge_{molA}_{molB}'
        dg_0 = json_0['estimate']['magnitude']
        file_1 = results_1 / file.split('/')[1]
        file_2 = results_2 / file.split('/')[1]
        with open(file_1) as stream:
            json_1 = json.load(stream)
        with open(file_2) as stream:
            json_2 = json.load(stream)
        dg_1 = json_1['estimate']['magnitude']
        dg_2 = json_2['estimate']['magnitude']
        dgs = [dg_0, dg_1, dg_2]
        dg = np.mean(dgs)
        std = np.std(dgs)
        if edge_name not in edges_dict:
            edges_dict[edge_name] = {
                    'ligand_a': molA,
                    'ligand_b': molB,
            }
        edges_dict[edge_name].update({runtype: [dg, std]})
    with open(output, 'w') as f:
        f.write("# Ligand1,Ligand2,calc_DDG,calc_dDDG(std)\n")
        for edge, data in edges_dict.items():
            ddg = data['complex'][0] - data['solvent'][0]
            error = np.sqrt(data['complex'][1]**2 + data['solvent'][1]**2)
            molA = data['ligand_a']
            molB = data['ligand_b']
            f.write(f"{molA},{molB},{ddg:.2f},{error:.2f}\n")


if __name__ == "__main__":
    extract()
