import pytest
import torch
from magnumnp import *
import numpy as np
import pathlib

import sys
sys.modules.pop("run", None)
import run

def test_slonczewski2():
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "log.dat"
    ref_path = this_dir / "ref" / "log2.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_x = torch.from_numpy(data[:, 4])
    data_y = torch.from_numpy(data[:, 5])
    data_z = torch.from_numpy(data[:, 6])

    ref_x = torch.from_numpy(ref[:, 4])
    ref_y = torch.from_numpy(ref[:, 5])
    ref_z = torch.from_numpy(ref[:, 6])

    torch.testing.assert_close(data_x, ref_x, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_y, ref_y, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_z, ref_z, atol=1e-3, rtol=1e-3)
