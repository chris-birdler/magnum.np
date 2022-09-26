import pytest
import pathlib
from magnumnp import *

#test CI
def test_read_vti():
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "fields" / "fields.vti"

    mesh, fields = read_vti(filename)
    assert len(fields) == 3
    assert fields['m'].shape == (100, 25, 1, 3)
    assert fields['h'].shape == (100, 25, 1, 3)
    assert fields['domain1'].shape == (100, 25, 1)

def test_write_vti():
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "fields" / "fields.vti"

    mesh, fields = read_vti(filename)
    assert len(fields) == 3
    assert fields['m'].shape == (100, 25, 1, 3)
    assert fields['h'].shape == (100, 25, 1, 3)
    assert fields['domain1'].shape == (100, 25, 1)
