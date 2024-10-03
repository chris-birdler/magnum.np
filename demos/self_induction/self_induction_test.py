import pytest
import torch
from magnumnp import *
from self_induction_run import run_self_induction
import numpy as np
import pathlib

def test_self_induction():
    run_self_induction()

    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "log.dat"
    ref_path = this_dir / "data" / "ref.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_x = torch.from_numpy(data[:, 1])
    data_y = torch.from_numpy(data[:, 2])
    data_z = torch.from_numpy(data[:, 3])

    ref_x = torch.from_numpy(ref[:, 1])
    ref_y = torch.from_numpy(ref[:, 2])
    ref_z = torch.from_numpy(ref[:, 3])

    torch.testing.assert_close(data_x, ref_x, atol=0, rtol=1e-2)
    torch.testing.assert_close(data_y, ref_y, atol=0, rtol=1e-2)
    torch.testing.assert_close(data_z, ref_z, atol=0, rtol=1e-2)

