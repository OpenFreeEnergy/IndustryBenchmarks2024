import numpy as np
import glob
import json
import gufe
from gufe.tokenization import JSON_HANDLER
from gufe import SmallMoleculeComponent as SMC
from cinnabar import Measurement, FEMap
from openff.units import unit
import pathlib
import click
import csv
from tqdm import tqdm


def get_names_from_unit_results(result) -> tuple[str, str]:
    # Result to tuple of ligand names
    # find string ligA to ligB repeat 0 generation 0
    # ligand names could have a space or an underscore in their name
    nm = list(result["unit_results"].values())[0]["name"]
    toks = nm.split(" to ")
    toks_2 = toks[1].split(" repeat")
    return toks[0], toks_2[0]


def get_names(result) -> tuple[str, str]:
    # get the name from the SmallMoleculeComponent
    list_of_pur = list(result["protocol_result"]["data"].values())[0]
    pur = list_of_pur[0]
    lig_A = pur["inputs"]["stateA"]["components"]["ligand"]
    lig_B = pur["inputs"]["stateB"]["components"]["ligand"]
    return SMC.from_dict(lig_A).name, SMC.from_dict(lig_B).name


def get_type(res):
    list_of_pur = list(res["protocol_result"]["data"].values())[0]
    pur = list_of_pur[0]
    components = pur["inputs"]["stateA"]["components"]

    if "solvent" not in components:
        return "vacuum"
    elif "protein" in components:
        return "complex"
    else:
        return "solvent"


def get_type_from_file_path(res):
    file_path = str(list(res["unit_results"].values())[0]["outputs"]["nc"])
    if "solvent" in file_path:
        return "solvent"
    elif "complex" in file_path:
        return "complex"
    # is anyone even running vacuum?
    elif "vacuum" in file_path:
        return "vaccum"
    else:
        raise ValueError("Can't guess simulation type")


def load_results(f):
    # path to deserialized results
    with open(f, "r") as fd:
        result = json.load(fd, cls=JSON_HANDLER.decoder)
    return result

def parse_ligand_network(
    input_ligand_network_path: str,
) -> gufe.LigandNetwork:
    with open(input_ligand_network_path) as f:
        graphml = f.read()

    network = gufe.LigandNetwork.from_graphml(graphml)
    return network


def get_dg(
        calc_data: dict[str, dict[str, float]],
        filename,
) -> None:
    """
    Helper method to write out MLE derived dG values.

    Parameters
    ----------
    calc_data: dict[str, dict[str, float]]
      The calculated DDG data.
    filename : pathlib.Path
      Pathlib object for saving the calculated DG data to a .tsv file.
    """
    fe_results = []

    for entry in calc_data:
        m = Measurement(
            labelA=calc_data[entry]['ligand_a'],
            labelB=calc_data[entry]['ligand_b'],
            DG=calc_data[entry]['ddG'] * unit.kilocalorie_per_mole,
            uncertainty=calc_data[entry][
                            'ddG_err'] * unit.kilocalorie_per_mole,
            computational=True
        )
        fe_results.append(m)

    # Feed into the FEMap object
    femap = FEMap()

    for entry in fe_results:
        femap.add_measurement(entry)

    femap.generate_absolute_values()
    df = femap.get_absolute_dataframe()
    df = df.iloc[:, :3]
    df.rename({'label': 'ligand'}, axis='columns', inplace=True)
    df.to_csv(filename, sep='\t', index=False, header=True)

    return


@click.command
@click.option(
    "--results_0",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("results_0"),
    required=True,
    help=(
        "Path to the directory that contains all result json files "
        "for repeat 0, default: results_0."
    ),
)
@click.option(
    "--results_1",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("results_1"),
    required=True,
    help=(
        "Path to the directory that contains all result json files "
        "for repeat 1, default: results_1."
    ),
)
@click.option(
    "--results_2",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("results_2"),
    required=True,
    help=(
        "Path to the directory that contains all result json files "
        "for repeat 2, default: results_2."
    ),
)
@click.option(
    "--input_ligand_network_file",
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./alchemicalNetwork/ligand_network.graphml"),
    required=True,
    help=(
        "Path to the used ligand network file."
        "(default: 'alchemicalNetwork/ligand_network.graphml')"
    ),
)
@click.option(
    "output",
    "-o",
    type=click.File(mode="w"),
    default=pathlib.Path("ddg.tsv"),
    required=True,
    help="Path to the output tsv file in which the DDG values are be stored. "
    "The output contains ligand names, DDG [kcal/mol] as the mean across "
    "three repeats and the standard deviation. Default: ddg.tsv.",
)
@click.option(
    'output_DG',
    '-o_DG',
    type=click.File(mode='w'),
    default=pathlib.Path('dg.tsv'),
    required=True,
    help="Path to the output tsv file in which the DG values are be stored. "
         "The output contains ligand names, the MLE derived DG [kcal/mol] "
         "values and associated uncertainties. Default: dg.tsv.",
)
def extract(results_0, results_1, results_2, input_ligand_network_file, output, output_DG):
    files_0 = glob.glob(f"{results_0}/*.json")
    files_1 = glob.glob(f"{results_1}/*.json")
    files_2 = glob.glob(f"{results_2}/*.json")

    list_of_files = [files_0, files_1, files_2]

    # Loading in the input ligand network to check if there are
    # edges that don't have a results .json file
    input_ligand_network = parse_ligand_network(input_ligand_network_file)

    click.echo("Checking files for errors and missing files")
    files_with_errors = []
    missing_files = []
    transform_dict = {"solvent": [], "complex": []}
    for i, file_list in tqdm(enumerate(list_of_files), desc="Looping over replicas"):
        solvent_list = []
        complex_list = []
        for file in tqdm(file_list, desc="Looping over transformations"):
            with open(file, "r") as fd:
                results = json.load(fd)

            # First we check if someone passed in a network_setup.json
            if results.get("__qualname__") == "AlchemicalNetwork":
                click.echo(f"{file} is a network_setup.json, removing from file list")
                file_list.remove(file)
                continue

            #  Now  we check if someome passed in an input json
            if results.get("__qualname__") == "Transformation":
                click.echo(f"{file} is an input json, skipping")
                file_list.remove(file)
                continue

            # Now we check if there are unit_results
            if "unit_results" not in results.keys():
                click.echo(f"{file} has no unit results")
                files_with_errors.append(file)
                continue

            # Now we check that we have a estimate and uncertainty
            # We print the traceback as well
            if results.get("estimate") is None or results.get("uncertainty") is None:
                click.echo(f"{file} has no estimate or uncertainty")
                proto_failures = [
                    k
                    for k in results["unit_results"].keys()
                    if k.startswith("ProtocolUnitFailure")
                ]
                for proto_failure in proto_failures:
                    click.echo("\n")
                    click.echo(results["unit_results"][proto_failure]["traceback"])
                    click.echo(results["unit_results"][proto_failure]["exception"])
                    click.echo("\n")
                files_with_errors.append(file)
                continue

            # Now create a list of all the Transformations present in the
            # result folders
            try:
                runtype = get_type(results)
            except KeyError:
                click.echo("Guessing run type from file path")
                runtype = get_type_from_file_path(results)
            try:
                molA, molB = get_names(results)
            except (KeyError, IndexError):
                try:
                    click.echo("Guessing names from simulation name")
                    molA, molB = get_names_from_unit_results(results)
                except (KeyError, IndexError):
                    raise ValueError("Failed to guess names")

            # Add ligand pairs of all transformations present to the dictionary
            if runtype == "solvent":
                solvent_list.append((molA, molB))
            if runtype == "complex":
                complex_list.append((molA, molB))

        # Check if there are transformations without results .json files
        for e in input_ligand_network.edges:
            molA = e.componentA.name
            molB = e.componentB.name
            if (molA, molB) not in solvent_list:
                missing_files.append(f"solvent_{molA}_{molB} repeat {i}")
            if (molA, molB) not in complex_list:
                missing_files.append(f"complex_{molA}_{molB} repeat {i}")

        # append to the global transforms dictionary
        transform_dict['solvent'].extend(solvent_list)
        transform_dict['complex'].extend(complex_list)

    # Before starting to raise errors, print out a summary of health checks
    click.echo("=" * 80)
    click.echo("Summary of checks for missing files and file with errors")
    click.echo("=" * 80)
    click.echo(
        "Total number of transformations that should have completed (both solvent and complex): "
        f"{len(input_ligand_network.edges)*2}"
    )

    if files_with_errors or missing_files:
        click.echo(
            "There are issues with some transformations, please contact"
            " the OpenFE team for next steps"
        )
        click.echo(
            "Total number of failed transformations: "
            f"{len(files_with_errors) + len(missing_files)}"
        )
        click.echo(
            "    - Number of transformations with errors in the result "
            f".json file: {len(files_with_errors)}"
        )
        click.echo(
            "    - Number of transformations without output "
            f"results files: {len(missing_files)}"
        )

        # Print details: Which files failed?
        if files_with_errors:
            click.echo("=" * 80)
            click.echo("Following result files contain errors:")
            for file in files_with_errors:
                click.echo(file)
            click.echo("=" * 80)

        if missing_files:
            click.echo("=" * 80)
            click.echo("There are no result files for the following transformations:")
            for file in missing_files:
                click.echo(file)
            click.echo("=" * 80)


        if files_with_errors:
            raise ValueError(
                "There are issues with these transformations, please contact the "
                "OpenFE team for next steps"
            )


        # if no errors, check for partial missing files
        unique_solvent_transforms = len(set(transform_dict["solvent"]))
        unique_complex_transforms = len(set(transform_dict["complex"]))

        if ((unique_solvent_transforms != unique_complex_transforms) or
            (unique_solvent_transforms != len(transform_dict["solvent"]) / 3) or
            (unique_complex_transforms != len(transform_dict["complex"]) / 3)):
            errmsg = (
                "Some files are partially missing for some edges. "
                "If you intentionally removed an edge, please check "
                "that you did not keep excess solvent or complex results JSON files."
            )
            raise ValueError(errmsg)


    # Check if there are .json files in the provided folders
    if len(files_0) == 0 or len(files_1) == 0 or len(files_2) == 0:
        errmsg = (
            "No .json files found in at least one of the results folders"
            ". Please check your specified file paths. Got folder names:"
            f" {results_0}, {results_1}, {results_2}."
        )
        raise ValueError(errmsg)

    else:
        click.echo("All simulations finished successfully!")
    click.echo("=" * 80)

    # Start extracting results
    click.echo("Start extracting results")
    edges_dict = dict()
    for file in files_0:
        click.echo(f"Reading file {file}")
        json_0 = load_results(file)
        try:
            runtype = get_type(json_0)
        except KeyError:
            click.echo("Guessing run type from file path")
            runtype = get_type_from_file_path(json_0)
        try:
            molA, molB = get_names(json_0)
        except (KeyError, IndexError):
            try:
                click.echo("Guessing names from simulation name")
                molA, molB = get_names_from_unit_results(json_0)
            except (KeyError, IndexError):
                raise ValueError("Failed to guess names")
        edge_name = f"edge_{molA}_{molB}"
        dg_0 = json_0["estimate"].m
        file_1 = results_1 / file.split("/")[-1]
        file_2 = results_2 / file.split("/")[-1]
        click.echo(f"Reading file {file_1}")
        json_1 = load_results(file_1)
        click.echo(f"Reading file {file_2}")
        json_2 = load_results(file_2)
        dg_1 = json_1["estimate"].m
        dg_2 = json_2["estimate"].m
        dgs = [dg_0, dg_1, dg_2]
        dg = np.mean(dgs)
        std = np.std(dgs)
        if edge_name not in edges_dict:
            edges_dict[edge_name] = {
                "ligand_a": molA,
                "ligand_b": molB,
            }
        edges_dict[edge_name].update({runtype: [dg, std]})

    writer = csv.writer(
        output,
        delimiter="\t",
        lineterminator="\n",  # to exactly reproduce previous, prefer "\r\n"
    )
    writer.writerow(["ligand_i", "ligand_j", "DDG(i->j) (kcal/mol)",
                     "uncertainty (kcal/mol)"])
    calc_data = {}
    for edge, data in edges_dict.items():
        ddg = data["complex"][0] - data["solvent"][0]
        error = np.sqrt(data["complex"][1] ** 2 + data["solvent"][1] ** 2)
        molA = data["ligand_a"]
        molB = data["ligand_b"]
        writer.writerow([molA, molB, round(ddg, 2), round(error, 2)])
        calc_data[edge] = {}
        calc_data[edge]['ligand_a'] = molA
        calc_data[edge]['ligand_b'] = molB
        calc_data[edge]['ddG'] = ddg
        calc_data[edge]['ddG_err'] = error

    get_dg(calc_data, output_DG)


if __name__ == "__main__":
    extract()
