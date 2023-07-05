import pytest
import torch
from magnumnp import *
from softmagnetic_composite_run import run_softmagnetic_composite
import numpy as np
import pathlib

def test_softmagnetic_composite():    
    run_softmagnetic_composite()
    
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "m.dat"
    ref_path = this_dir / "ref" / "m_test.dat"
    
    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)
    
    data_h = torch.from_numpy(data[:,1])
    data_m = torch.from_numpy(data[:,4])
    
    ref_h = torch.from_numpy(ref[:,1])
    ref_m = torch.from_numpy(ref[:,4])
    
    torch.testing.assert_close(data_h, ref_h, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_m, ref_m, atol=1e-3, rtol=1e-3)
