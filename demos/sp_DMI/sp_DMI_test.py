import pytest
import torch
from magnumnp import *
from sp_DMI_run import run_sp_DMI
import numpy as np
import pathlib

def test_sp_DMI():
    run_sp_DMI()

    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "m0_magnumnp.dat"
    ref_path = this_dir / "ref" / "m_test.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_x = torch.from_numpy(data[1:-1,1])
    data_z = torch.from_numpy(data[1:-1:2,3])
    ref_x = torch.from_numpy(ref[1:-1,1])
    ref_z = torch.from_numpy(ref[1:-1:2,3])

    torch.testing.assert_close(data_x, ref_x, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_z, ref_z, atol=1e-3, rtol=1e-3)
