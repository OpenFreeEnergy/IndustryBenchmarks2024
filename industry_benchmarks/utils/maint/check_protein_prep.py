import numpy as np
import pathlib
import click
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


@click.command
@click.option(
    '--prepared',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared PDB file",
)
@click.option(
    '--original',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the original PDB file",
)
def run(prepared, original):
    """
    Compare the prepared PDB with its original form, checking that there
    are not major differences in atom counts and positions.

    Parameters
    ----------
    prepared : pathlib.Path
      A Path to the prepared PDB file.
    original : pathlib.Path
      A Path to the `original` PDB file.
    """
    prep = mda.Universe(prepared)
    prev = mda.Universe(original)

    prep_prot = prep.select_atoms('protein and not resname ACE NME NMA')
    prev_prot = prev.select_atoms('protein and not resname ACE NME NMA')

    for resi, resj in zip(prep_prot.residues, prev_prot.residues):
        if len(resi.atoms) != len(resj.atoms):
            print("number of atoms does not match: ", resi.resname, resi.resid)
        else:
            arr = np.abs(resi.atoms.positions - resj.atoms.positions)
            if np.any(arr > 1.0):
                print(">1 A deviation in residue: ", resi.resname, resi.resid)

    prep_nonprot = prep.select_atoms('not protein and not resname NMA')
    prev_nonprot = prev.select_atoms('not protein and not resname NMA')

    # Check that the number of atoms match for non-protein residues
    if len(prep_nonprot.atoms) != len(prev_nonprot.atoms):
        print("number of non protein atoms does not match")

    # Inspect disulfide bridges
    # first guess bonds
    prep.atoms.guess_bonds()
    sgs = prep.select_atoms('name SG and resname CYS CYX')
    distances = distance_array(sgs.atoms.positions, sgs.atoms.positions)

    for i in range(len(sgs)):
        for j in range(len(sgs)):
            # skip self
            if i == j:
                continue
            else:
                if distances[i, j] < 2.5:
                    if sgs.atoms[j] not in sgs.atoms[i].bonded_atoms:
                        if 'H' in sgs.atoms[i].bonded_atoms.elements:
                            print(f"ERROR: possible missed bridge between {sgs.atoms[j]} and {sgs.atoms[i]}")
                        else:
                            print(f"WARNING: unbonded disulfide bridge found: {sgs.atoms[j]} and {sgs.atoms[i]}")
                    else:
                        print(f"LOG: disulfide bridge found {sgs.atoms[j]} and {sgs.atoms[i]}")


if __name__ == "__main__":
    run()
