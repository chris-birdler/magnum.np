import pytest
import torch
from magnumnp import *
import numpy as np
import pathlib

import sys
sys.modules.pop("run", None)
import run

def test_rkky():    
    Hk = 2 * 1e5 / 1.
    
    this_dir = pathlib.Path(__file__).resolve().parent
    data_path = this_dir / "data" / "log.dat"
    ref_path = this_dir / "ref" / "log_test.dat"

    data = np.loadtxt(data_path)
    ref = np.loadtxt(ref_path)

    data_m = torch.from_numpy(data[:, 3])
    data_h = torch.from_numpy(data[:, 6]/Hk)
    
    ref_m = torch.from_numpy(ref[:, 3])
    ref_h = torch.from_numpy(ref[:, 6]/Hk)
    
    torch.testing.assert_close(data_m, ref_m, atol=1e-3, rtol=1e-3)
    torch.testing.assert_close(data_h, ref_h, atol=1e-3, rtol=1e-3) 
