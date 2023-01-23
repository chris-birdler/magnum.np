magnum.np
=========

[![Documentation Status](https://readthedocs.com/projects/magnumnp-magnumnp/badge/?version=latest&token=14bfc008bc74536d868c875c61777aea0733cb782c08690a29aecf5761e21f1f)](https://magnumnp-magnumnp.readthedocs-hosted.com/en/latest/?badge=latest)

magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference method. It implements state-of-the-art algorithms and is based on [pytorch](http://www.pytorch.org/), which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in Python which leads to very readable yet flexible code. Due to [pytorch](http://www.pytorch.org/) integration, extensive postprocessing can be done directly in the simulations scripts. Furthermore [pytorch](http://www.paraview.org/)'s autograd feature makes it possible to solve inverse problems without significant modifications of the code. This manual is meant to give you both a quick start and a reference to magnum.np.


Features
--------
* Explicit / Implicit time-integration of the Landau-Lifshitz-Gilbert Equation
* Fast FFT Demagnetization-field computation optimized for small memory footprint
* Fast FFT Oersted-field optimized for small memory footprint
* Arbitrary Material Parameters variing in space and time
* Spin-torque model by Zhang and Li, Slonczewski
* Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface, bulk, D2d)
* String method for energy barrier computations
* Sophisticated domain handling, e.g. for spatially varying material parameters
* Seemingless VTK import / export via [pyvista](https://docs.pyvista.org/)
* Inverse Problems via [pytorch](http://www.paraview.org/)'s autograd feature


Example
-------
The following demo code shows the solution of the MuMag Standard Problem #5 and can be found in the demos directory:

```python
from magnumnp import *
import torch

Timer.enable()

# initialize state
n  = (40, 40, 1)
dx = (2.5e-9, 2.5e-9, 10e-9)
mesh = Mesh(n, dx)

state = State(mesh)
state.material = {
    "Ms": 8e5,
    "A": 1.3e-11,
    "alpha": 0.1,
    "xi": 0.05,
    "b": 72.17e-12
    }

# initialize magnetization that relaxes into s-state
state.m = state.Constant([0,0,0])
state.m[:20,:,:,1] = -1.
state.m[20:,:,:,1] = 1.
state.m[20,20,:,1] = 0.
state.m[20,20,:,2] = 1.

state.j = state.Tensor([1e12, 0, 0])

# initialize field terms
demag    = DemagField()
exchange = ExchangeField()
torque   = SpinTorqueZhangLi()

# initialize sstate
llg = LLGSolver([demag, exchange])
llg.relax(state)
write_vti(state.m, "data/m0.vti", state)

# perform integration with spin torque
llg = LLGSolver([demag, exchange, torque])
logger = ScalarLogger("data/m.dat", ['t', 'm'])
while state.t < 5e-9:
    llg.step(state, 1e-10)
    logger << state

Timer.print_report()
```

Documentation
-------------
The documentation is located in the doc directory and can be built using [sphinx](https://www.sphinx-doc.org).
For example the following commands build an HTML documentation of the actual source code.
```
cd doc
make html
```

Alternatively, the latest version of the documentation is always available on [https://magnumnp-magnumnp.readthedocs-hosted.com/en/latest/](https://magnumnp-magnumnp.readthedocs-hosted.com/en/latest/)


Citation
--------
If you use magnum.np in your work or publication, please cite the following reference:

[1] Bruckner, Florian, et al. "magnum.np -- A pytorch based GPU enhanced Finite Difference Micromagnetic Simulation Framework for High Level Development and Inverse Design", to be published (2023).


Contributing
------------
Contributions are gratefully accepted.
The source code is hosted on [www.gitlab.com/magnum.np/magnum.np](www.gitlab.com/magnum.np/magnum.np).
If you have any issues or question, just open an issue via gitlab.com.
To contribute code, fork our repository on gitlab.com and create a corresponding merge request.
