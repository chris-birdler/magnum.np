import pytest
import torch
from magnumnp import *
from sp_domainwall_pinning_run import run_sp_domainwall_pinning
import numpy as np
import pathlib

def test_domainwall_pinning():
    run_sp_domainwall_pinning()

    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir /"data"/"m.dat"
    ref_path = this_dir /"ref"/"m_test.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_h = torch.from_numpy(data[:,2])
    data_m = torch.from_numpy(data[:,5])

    ref_h = torch.from_numpy(ref[:,2])
    ref_m = torch.from_numpy(ref[:,5])

    torch.testing.assert_close(data_h, ref_h, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_m, ref_m, atol=1e-3, rtol=1e-3)
