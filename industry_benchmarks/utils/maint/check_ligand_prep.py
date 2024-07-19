from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pathlib
import click


@click.command
@click.option(
    '--ligands',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared SDF file",
)
def run(ligands):
    """
    Inspect the prepared ligands file, making sure that there are no
    duplicates.

    Parameters
    ----------
    ligands : pathlib.Path
      A path to the prepared SDF file.
    """
    uncharger = rdMolStandardize.Uncharger()
    supp = Chem.SDMolSupplier(ligands)
    umols = [uncharger.uncharge(m) for m in supp]
    inchis = [Chem.MolToInchi(m) for m in umols]
    
    print("number of ligands: ", len(inchis))

    duplicates = []

    for key1, i in enumerate(inchis):
        for key2, j in enumerate(inchis):
            if key1 == key2:
                continue
            if i == j:
                duplicates.append(j)
                inchis.pop(key2)


    if len(duplicates) > 0:
        print("Error: Number of duplicate ligands found: ", len(duplicates))
    else:
        print("Info: No duplicate ligands found")


if __name__ == "__main__":
    run()
