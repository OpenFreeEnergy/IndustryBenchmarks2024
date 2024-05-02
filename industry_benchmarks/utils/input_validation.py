import string
import click
import tempfile
import pathlib
import gufe
from openff.units import unit
import openfe
from openfe.protocols.openmm_md.plain_md_methods import PlainMDProtocol
from rdkit import Chem


def get_settings():
    """
    Utility method for getting MDProtocol settings.

    These settings mostly follow defaults but use very short
    simulation times to avoid being too much of a burden on users' machines.
    """
    settings = openfe.protocols.openmm_md.plain_md_methods.PlainMDProtocol.default_settings()
    settings.simulation_settings.equilibration_length_nvt = 1 * unit.picosecond
    settings.simulation_settings.equilibration_length = 1 * unit.picosecond
    settings.simulation_settings.production_length = 1 * unit.picosecond
    return settings


def run_md(dag):
    """
    Run a DAG and check it was ok.

    Parameters
    ----------
    dag : openfe.ProtocolDAG
      A ProtocolDAG to execute.

    Raises
    ------
    AssertionError
      If any of the simulation Units failed.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = pathlib.Path(tmpdir)
        dagres = gufe.protocols.execute_DAG(
            dag,
            shared_basedir=workdir,
            scratch_basedir=workdir,
            keep_shared=False,
            raise_error=True,
            n_retries=0,
        )

    # If everything is good then tell the user
    assert dagres.ok()
    print("SIMULATION COMPLETE")


@click.command
@click.option(
    '--pdb',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared PDB file to validate",
)
@click.option(
    '--cofactors',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=None,
    help="Path to the prepared cofactors SDF file (optional)",
)
def run_inputs(pdb, cofactors):
    """
    Validate input files by running a short MD simulation

    Parameters
    ----------
    pdb : pathlib.Path
      A Path to a protein PDB file.
    cofactors : Optional[pathlib.Path]
      A Path to an SDF file containing the system's cofactors.
    """
    # Create the solvent and protein components
    solv = openfe.SolventComponent()
    prot = openfe.ProteinComponent.from_pdb_file(str(pdb))

    # Store there in a components dictionary
    components_dict = {
        'protein': prot,
        'solvent': solv,
    }

    # If we have cofactors, populate them and store them based on
    # an single letter index (we assume no more than len(alphabet) cofactors)
    if cofactors is not None:
        cofactors = [
            openfe.SmallMoleculeComponent(m)
            for m in Chem.SDMolSupplier(str(cofactors), removeHs=False)
        ]

        for cofactor, entry in zip(cofactors, string.ascii_lowercase):
            components_dict[entry] = cofactor

    # Create the ChemicalSystem
    system = openfe.ChemicalSystem(components_dict)

    # Get the settings and create the protocol
    settings = get_settings()
    protocol = PlainMDProtocol(settings=settings)

    # Now create the DAG and run it
    dag = protocol.create(stateA=system, stateB=system, mapping=None)
    run_md(dag)


if __name__ == "__main__":
    run_inputs()
