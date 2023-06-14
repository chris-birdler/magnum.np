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

run_FMR = run_module.run_FMR

def test_FMR():    
    
    run_FMR()
    
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir /"data"/"log.dat"
    ref_path = this_dir /"ref"/"log_test.dat"
    
    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)
    
    data_y = torch.from_numpy(data[:, 2])
    ref_y = torch.from_numpy(ref[:, 2])
    
    torch.testing.assert_close(data_y, ref_y, atol=1e-10, rtol=1e-10) 