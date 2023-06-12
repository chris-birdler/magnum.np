import pytest
import torch
from magnumnp import *
from .run import run_sp4
import numpy as np
    
def test_sp4():    
    
    run_sp4()
    
    data = np.loadtxt("sp4/data/log.dat")
    ref = np.loadtxt("sp4/ref/m_test.dat")
    
    data_x = torch.from_numpy(data[:, 1])
    data_y = torch.from_numpy(data[:, 2])
    data_z = torch.from_numpy(data[:, 3])
    
    ref_x = torch.from_numpy(ref[:, 1])
    ref_y = torch.from_numpy(ref[:, 2])
    ref_z = torch.from_numpy(ref[:, 3])
    
    torch.testing.assert_close(data_x, ref_x, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(data_y, ref_y, atol=1e-10, rtol=1e-10) 
    torch.testing.assert_close(data_z, ref_z, atol=1e-10, rtol=1e-10) 
    
    