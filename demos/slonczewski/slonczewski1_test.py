import pytest
import torch
from magnumnp import *
from slonczewski1_run import run_slonczewski1
import numpy as np
import pathlib

def test_slonczewski1():
    run_slonczewski1()

    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data1" / "log.dat"
    ref_path = this_dir / "ref" / "log1.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_x = torch.from_numpy(data[:, 1])
    data_y = torch.from_numpy(data[:, 2])
    data_z = torch.from_numpy(data[:, 3])

    ref_x = torch.from_numpy(ref[:, 1])
    ref_y = torch.from_numpy(ref[:, 2])
    ref_z = torch.from_numpy(ref[:, 3])

    torch.testing.assert_close(data_x, ref_x, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_y, ref_y, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_z, ref_z, atol=1e-3, rtol=1e-3)
