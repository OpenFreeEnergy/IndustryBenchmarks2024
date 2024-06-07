import pathlib
from openmmtools import multistate
import numpy as np
import MDAnalysis as mda
from openfe_analysis import FEReader


def compute_mbar_energies(analyzer):
    """
    Returns
    ------
    u_ln: energy matrix of shape (K,N) indexed by k,n
        K is the total number of states observables are desired.
        N is the total number of samples drawn from ALL states.
        The nth configuration is the energy evaluated in the kth thermodynamic state.
    N_l: 1-D iterable of shape K
        The number of samples drawn from each kth state.
    """
    analyzer.use_full_trajectory = True
    u_ln, N_l = analyzer._compute_mbar_decorrelated_energies()
    return u_ln, N_l


def get_replica_state_indices(analyzer):
    energy_data = list(analyzer._read_energies(truncate_max_n_iterations=True))
    replicas_state_indices = energy_data[-1]
    return replicas_state_indices


def subsample_traj(simulation, hybrid_system_pdb, lambda_windows, outfile):
    for i in range(0, lambda_windows):
        u = mda.Universe(hybrid_system_pdb, simulation,
                         format=FEReader, state_id=i)
        skip = round(len(u.trajectory) / 20)
        out_traj = pathlib.Path(f'{outfile}_{i}.xtc')
        with mda.Writer(str(out_traj), n_atoms=len(u.atoms)) as w:
            for ts in u.trajectory[::skip]:
                w.write(u.atoms)


def extract_data(simulation, checkpoint, hybrid_pdb, outfile, out_traj='out'):
    """
    Extract the MBAR-ready energy matrix and replica state indices.

    Parameters
    ----------
    simulation: pathlib.Path
        Path to the simulation `.nc` file
    checkpoint: pathlib.Path
        Path to the checkpoint `.chk` file
    hybrid_pdb: pathlib.Path
        Path to the `.pdb` file of the hybrid system
    outfile: pathlib.Path
        Path to the output `.npz` file
    out_traj: pathlib.Path
        Path to the output `.xtc` files. A separate file is created for every
        lambda window. The state number is appended to the filename. Default: 'out'
    """
    if not simulation.is_file() or not checkpoint.is_file():
        errmsg = "Either the simulation or checkpoint file could not be found"
        raise ValueError(errmsg)

    reporter = multistate.MultiStateReporter(
        storage=simulation.as_posix(),
        open_mode='r',
        checkpoint_storage=checkpoint.as_posix()
    )

    analyzer = multistate.MultiStateSamplerAnalyzer(reporter)
    u_ln, N_l = compute_mbar_energies(analyzer)
    replicas_state_indices = get_replica_state_indices(analyzer)
    np.savez(outfile, u_ln=u_ln, N_l=N_l, replicas_state_indices=replicas_state_indices)

    lambda_windows = len(replicas_state_indices)
    subsample_traj(simulation, hybrid_pdb, lambda_windows, out_traj)
