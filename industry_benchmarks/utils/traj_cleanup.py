import pathlib
from openmmtools import multistate
import numpy as np


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


def extract_data(simulation, checkpoint, outfile):
    """
    Extract the MBAR-ready energy matrix and replica state indices.

    Parameters
    ----------
    simulation: pathlib.Path
        Path to the simulation `.nc` file
    checkpoint: pathlib.Path
        Path to the checkpoint `.chk` file
    outfile: pathlib.Path
        Path to the output `.npz` file
    """

    reporter = multistate.MultiStateReporter(
        storage=simulation.as_posix(),
        open_mode='r',
        checkpoint_storage=checkpoint.as_posix()
    )
    analyzer = multistate.MultiStateSamplerAnalyzer(reporter)
    u_ln, N_l = compute_mbar_energies(analyzer)
    replicas_state_indices = get_replica_state_indices(analyzer)
    np.savez(outfile, u_ln=u_ln, N_l=N_l, replicas_state_indices=replicas_state_indices)
