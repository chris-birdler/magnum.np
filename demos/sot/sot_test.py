import pytest
import torch
from magnumnp import *
from .run import run_sot
import numpy as np
    
def test_sot():    
    
    run_sot()
    
    data = np.loadtxt("sot/data/log.dat")
    ref = np.loadtxt("sot/ref/ref_test.dat")
    
    data_x = torch.from_numpy(data[:, 1])
    data_y = torch.from_numpy(data[:, 2])
    
    ref_x = torch.from_numpy(ref[:, 1])
    ref_y = torch.from_numpy(ref[:, 2])
    
    torch.testing.assert_close(data_x, ref_x, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(data_y, ref_y, atol=1e-10, rtol=1e-10) 
