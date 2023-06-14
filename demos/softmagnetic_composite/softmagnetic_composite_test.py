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

run_softmagnetic_composite = run_module.run_softmagnetic_composite
    
def test_softmagnetic_composite():    
    
    run_softmagnetic_composite()
    
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir /"data"/"m.dat"
    ref_path = this_dir /"ref"/"m_test.dat"
    
    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)
    
    data_h = torch.from_numpy(data[:, 1])
    data_m = torch.from_numpy(data[:, 4])
    
    ref_h = torch.from_numpy(ref[:, 1])
    ref_m = torch.from_numpy(ref[:, 4])
    
    torch.testing.assert_close(data_h, ref_h, atol=1e-10, rtol=1e-10)
    torch.testing.assert_close(data_m, ref_m, atol=1e-10, rtol=1e-10)
    