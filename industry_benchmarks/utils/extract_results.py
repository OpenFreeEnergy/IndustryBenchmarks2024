import numpy as np
import glob
import json

files_0 = glob.glob("results_0/*.json")
files_1 = glob.glob("results_1/*.json")
files_2 = glob.glob("results_2/*.json")
if len(files_0) != len(files_1) != len(files_2):
    print(f'CAVE: Some edges did not finish: repeat 0: {len(files_0)} edges, '
          f'repeat 1: {len(files_1)} edges, repeat 2: {len(files_2)} edges')


def get_names(result) -> tuple[str, str]:
    # Result to tuple of ligand names
    nm = list(result['unit_results'].values())[0]['name']
    toks = nm.split()
    if toks[2] == 'repeat':
        return toks[0], toks[1]
    else:
        return toks[0], toks[2]


edges_dict = dict()
for file in files_0:
    with open(file) as stream:
        results = json.load(stream)
    runtype = file.split('/')[1].split('_')[-2]
    molA, molB = get_names(results)
    edge_name = f'edge_{molA}_{molB}'
    dg_1 = results['estimate']['magnitude']
    file_2 = 'results_1/%s'%file.split('/')[1]
    file_3 = 'results_2/%s'%file.split('/')[1]
    with open(file_2) as stream:
        results_2 = json.load(stream)
    with open(file_3) as stream:
        results_3 = json.load(stream)
    dg_2 = results_2['estimate']['magnitude']
    dg_3 = results_3['estimate']['magnitude']
    dgs = [dg_1, dg_2, dg_3]
    dg = np.mean(dgs)
    std = np.std(dgs)
    if edge_name not in edges_dict:
        edges_dict[edge_name] = {
                'ligand_a': molA,
                'ligand_b': molB,
        }
    edges_dict[edge_name].update({runtype: [dg, std]})
with open('ddg.csv', 'w') as f:
    f.write("# Ligand1,Ligand2,calc_DDG,calc_dDDG(std)\n")
    for edge, data in edges_dict.items():
        ddg = data['complex'][0] - data['solvent'][0]
        error = np.sqrt(data['complex'][1]**2 + data['solvent'][1]**2)
        molA = data['ligand_a']
        molB = data['ligand_b']
        f.write(f"{molA},{molB},{ddg:.2f},{error:.2f}\n")
