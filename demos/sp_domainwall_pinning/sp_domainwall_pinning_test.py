import pytest
import torch
from magnumnp import *
from .run import run_domainwall_pinning
import numpy as np
    
def test_domainwall_pinning():    
    
    run_domainwall_pinning()

    data = np.loadtxt("sp_domainwall_pinning/data/m.dat")
    ref = np.loadtxt("sp_domainwall_pinning/ref/m_test.dat")
    
    data_h = torch.from_numpy(data[:,2])
    data_m = torch.from_numpy(data[:,5])
    ref_h = torch.from_numpy(ref[:,2])
    ref_m = torch.from_numpy(ref[:,5])
    
    torch.testing.assert_close(data_h, ref_h, atol=1e-10, rtol=1e-10) 
    torch.testing.assert_close(data_m, ref_m, atol=1e-10, rtol=1e-10)