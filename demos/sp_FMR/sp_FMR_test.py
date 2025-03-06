import pytest
import torch
from magnumnp import *
import run
import numpy as np
import pathlib

def test_FMR():
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "log.dat"
    ref_path = this_dir / "ref" / "log_test.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_y = torch.from_numpy(data[:, 2])
    ref_y = torch.from_numpy(ref[:, 2])

    torch.testing.assert_close(data_y, ref_y, atol=1e-3, rtol=1e-3)
