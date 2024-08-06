import gufe
import numpy as np
import glob
import json
from gufe.tokenization import JSON_HANDLER
from gufe import SmallMoleculeComponent as SMC
import pathlib
import click
import csv


def get_names(result) -> tuple[str, str]:
    # Result to tuple of ligand names
    # find string ligA to ligB repeat 0 generation 0
    # ligand names could have a space or an underscore in their name
    nm = list(result['unit_results'].values())[0]['name']
    toks = nm.split(' to ')
    toks_2 = toks[1].split(' repeat')
    return toks[0], toks_2[0]

def get_names(result) -> tuple[str, str]:
    # get the name from the SmallMoleculeComponent
    list_of_pur = list(result['protocol_result']['data'].values())[0]
    pur = list_of_pur[0]
    lig_A = pur['inputs']['stateA']['components']['ligand']
    lig_B = pur['inputs']['stateB']['components']['ligand']
    return SMC.from_dict(lig_A).name, SMC.from_dict(lig_B).name

def get_type(res):
    list_of_pur = list(res['protocol_result']['data'].values())[0]
    pur = list_of_pur[0]
    components = pur['inputs']['stateA']['components']

    if 'solvent' not in components:
        return 'vacuum'
    elif 'protein' in components:
        return 'complex'
    else:
        return 'solvent'

def load_results(f):
    # path to deserialized results
    result = json.load(open(f, 'r'), cls=JSON_HANDLER.decoder)
    if result['estimate'] is None or result['uncertainty'] is None:
        click.echo(f"Calculations for {f} did not finish successfully!")
        proto_failures = [k for k in result["unit_results"].keys() if k.startswith("ProtocolUnitFailure")]
        for proto_failure in proto_failures:
            click.echo("\n")
            click.echo(results["unit_results"][proto_failure]["traceback"])
            click.echo(results["unit_results"][proto_failure]["exception"])
            click.echo("\n")
        raise ValueError("Calculations did not finish successfully")

    return result


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
    'output',
    '-o',
    type=click.File(mode='w'),
    default=pathlib.Path('ddg.tsv'),
    required=True,
    help="Path to the output tsv file in which the DDG values are be stored. "
         "The output contains ligand names, DDG [kcal/mol] as the mean across "
         "three repeats and the standard deviation. Default: ddg.tsv.",
)
def extract(results_0, results_1, results_2, output):
    files_0 = glob.glob(f"{results_0}/*.json")
    files_1 = glob.glob(f"{results_1}/*.json")
    files_2 = glob.glob(f"{results_2}/*.json")
    # Check if there are .json files in the provided folders
    if len(files_0) == 0 or len(files_1) == 0 or len(files_2) == 0:
        errmsg = ('No .json files found in at least one of the results folders'
                  '. Please check your specified file paths. Got folder names:'
                  f' {results_0}, {results_1}, {results_2}.')
        raise ValueError(errmsg)
    # Check if there are missing files
    all_jsons = ([x.split('/')[-1] for x in files_0]
                 + [x.split('/')[-1] for x in files_1]
                 + [x.split('/')[-1] for x in files_2])
    missing_files = [x for x in set(all_jsons) if all_jsons.count(x) <= 2]
    if len(missing_files) > 0:
        errmsg = ('Some calculations did not finish and did not output a '
                  'result .json file. Number of .json files for the three '
                  f'repeats are: repeat 0: {len(files_0)} files, repeat 1: '
                  f'{len(files_1)} files, and repeat 2: {len(files_2)} files. '
                  f'Missing results have been found for {missing_files}.')
        raise ValueError(errmsg)

    # Start extracting results
    edges_dict = dict()
    for file in files_0:
        click.echo(f'Reading file {file}')
        json_0 = load_results(file)
        runtype = get_type(json_0)
        try:
            molA, molB = get_names(json_0)
        except (KeyError, IndexError):
            raise ValueError("Failed to guess names")
        edge_name = f'edge_{molA}_{molB}'
        dg_0 = json_0['estimate'].m
        file_1 = results_1 / file.split('/')[-1]
        file_2 = results_2 / file.split('/')[-1]
        click.echo(f'Reading file {file_1}')
        json_1 = load_results(file_1)
        click.echo(f'Reading file {file_2}')
        json_2 = load_results(file_2)
        dg_1 = json_1['estimate'].m
        dg_2 = json_2['estimate'].m
        dgs = [dg_0, dg_1, dg_2]
        dg = np.mean(dgs)
        std = np.std(dgs)
        if edge_name not in edges_dict:
            edges_dict[edge_name] = {
                    'ligand_a': molA,
                    'ligand_b': molB,
            }
        edges_dict[edge_name].update({runtype: [dg, std]})

    writer = csv.writer(
        output,
        delimiter="\t",
        lineterminator="\n",  # to exactly reproduce previous, prefer "\r\n"
    )
    writer.writerow(["ligand_i", "ligand_j", "DDG(i->j) (kcal/mol)",
                     "uncertainty (kcal/mol)"])
    for edge, data in edges_dict.items():
        ddg = data['complex'][0] - data['solvent'][0]
        error = np.sqrt(data['complex'][1] ** 2 + data['solvent'][1] ** 2)
        molA = data['ligand_a']
        molB = data['ligand_b']
        writer.writerow([molA, molB, round(ddg, 2), round(error, 2)])


if __name__ == "__main__":
    extract()
