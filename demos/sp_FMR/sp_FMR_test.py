import pytest
import torch
from magnumnp import *
from .run import run_FMR
import numpy as np
    
def test_FMR():    
    
    run_FMR()

    data = np.loadtxt("sp_FMR/data/log.dat")
    ref = np.loadtxt("sp_FMR/ref/m_test.dat")
    
    data_y = torch.from_numpy(data[:, 2])
    ref_y = torch.from_numpy(ref[:, 2])
    
    torch.testing.assert_close(data_y, ref_y, atol=1e-10, rtol=1e-10) 