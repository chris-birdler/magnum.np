import pytest
import torch
from magnumnp import *
import sys
import os
import numpy as np
import pathlib

current_dir = os.path.dirname(os.path.abspath(__file__))

parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

sys.path.insert(0, parent_dir)

subdirectory = os.path.basename(current_dir)

module_name = "{}.run".format(subdirectory)

run_module = __import__(module_name, fromlist=['run'])

run_DMI = run_module.run_DMI
    
def test_DMI():    
    
    run_DMI()
    
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir /"data"/"m0_magnumnp.dat"
    ref_path = this_dir /"ref"/"m_test.dat"
    
    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)
    
    data_x = torch.from_numpy(data[1:-1,1])
    data_z = torch.from_numpy(data[1:-1:2,3])
    ref_x = torch.from_numpy(ref[1:-1,1])
    ref_z = torch.from_numpy(ref[1:-1:2,3])
    
    torch.testing.assert_close(data_x, ref_x, atol=1e-6, rtol=0.02) 
    torch.testing.assert_close(data_z, ref_z, atol=1e-6, rtol=0.02)