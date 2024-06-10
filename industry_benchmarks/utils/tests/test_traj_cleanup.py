import pathlib

import MDAnalysis as mda
import pytest
from importlib import resources
from ..traj_cleanup import extract_data
import numpy as np


@pytest.fixture
def simulation():
    with resources.files('utils.tests.data.example_traj') as d:
        yield d / 'simulation.nc'


@pytest.fixture
def checkpoint():
    with resources.files('utils.tests.data.example_traj') as d:
        yield d / 'checkpoint.chk'


@pytest.fixture
def pdb():
    with resources.files('utils.tests.data.example_traj') as d:
        yield d / 'hybrid_system.pdb'


@pytest.fixture
def outfile():
    with resources.files('utils.tests.data.example_traj') as d:
        yield d / 'data.npz'


def test_extract_data(simulation, checkpoint, pdb, outfile, tmp_path):
    out_traj = tmp_path / 'out'
    out_npz = tmp_path / 'test.npz'

    extract_data(simulation, checkpoint, pdb, out_npz, out_traj)
    # Check that a .npz file was written out
    assert out_npz.is_file()
    # Check that the output is the same as the stored data
    test_output = np.load(tmp_path / 'test.npz')
    stored_output = np.load(outfile)
    assert list(stored_output) == list(test_output)
    u_ln = test_output['u_ln']
    N_l = test_output['N_l']
    replicas_state_indices = test_output['replicas_state_indices']
    # u_ln array should have a length of 11 (number of lambda windows)
    assert len(u_ln) == 11
    # N_l array should have a length of 11 (number of lambda windows)
    assert len(N_l) == 11
    # All values in N_l should be equal since we didn't subsample the data
    assert list(N_l).count(N_l[0]) == len(N_l)
    # Check the length of the u_ln matrix
    assert [len(n) == len(u_ln) * N_l[0] for n in u_ln]
    # Check that replicate state indices are of length N*L
    assert len(replicas_state_indices) == 11
    assert [len(r) == N_l[0] for r in replicas_state_indices]

    # Check output trajectories
    for i in range(11):
        # Check that trajectory files exist
        traj = pathlib.Path(f'{out_traj}_{i}.xtc')
        assert traj.is_file()
        u = mda.Universe(pdb, traj)
        # Check that we saved 20 frames
        assert len(u.trajectory) == 20
