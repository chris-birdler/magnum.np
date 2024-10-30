![magnum.np](logo.png)

magnum.np 2.0.1
===============
magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference method. It implements state-of-the-art algorithms and is based on [pytorch](http://www.pytorch.org/), which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in Python which leads to very readable yet flexible code. Due to [pytorch](http://www.pytorch.org/) integration, extensive postprocessing can be done directly in the simulations scripts. Furthermore [pytorch](http://www.paraview.org/)'s autograd feature makes it possible to solve inverse problems without significant modifications of the code. This manual is meant to give you both a quick start and a reference to magnum.np.

**Version 2.0:** The magnum.np interface slightly changed since version 2.0. Find more details about the necessary changes and the motivation [here](docs/changes.rst).

Features
--------
* Explicit / Implicit time-integration of the Landau-Lifshitz-Gilbert Equation
* Fast FFT Demagnetization-field computation optimized for small memory footprint
* Fast FFT Oersted-field optimized for small memory footprint
* Fast FFT Vector-Potential optimized for small memory footprint
* Periodic Boundary Conditions in 1D, 2D, and 3D (True and Pseudo-Periodic)
* Non-Equidistant Mesh for Multilayer Structures
* Arbitrary Material Parameters varying in space and time
* Spin-torque model by Slonczewski
* Spin-torque model by Zhang and Li
* Spin-Orbit torque (SOT)
* Bilinear and biquadratic Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface, bulk, D2d)
* String method for energy barrier computations
* [Eigenmode Solver](docs/eigensolver.rst) for efficient calculation of normal modes
* Sophisticated domain handling, e.g. for spatially varying material parameters
* efficient [Voronoi Code](docs/voronoi.rst) for 2D and 3D problems (including intergrain phase)
* Seamless VTK import / export via [pyvista](https://docs.pyvista.org/)
* Inverse Problems via [pytorch](http://www.pytorch.org/)'s autograd feature


Documented Demos:
-----------------
Demo scripts for various applications are available in the [demo](demos/README.md) directory.

The following demos are also stored on Google Colab, where they can directly be run without any local installation:

   * [Slonczewski Spin Torque](demos/slonczewski/run2.ipynb) ([Colab](https://colab.research.google.com/drive/1KJCd-mnTaveZruLV37OWzYhtuDGuu3rH))
   * [Softmagnetic Composite](demos/softmagnetic_composite/run.ipynb) ([Colab](https://colab.research.google.com/drive/1HazB7ydSYZKbtrQoPc9xE3U0d7uc-1Ir))
   * [Spin Orbit Torque](demos/sot/run.ipynb) ([Colab](https://colab.research.google.com/drive/1OWMH0_qqxM73rB5gK5pi7nFRtO4nO_N8))
   * [Standard Problem #4](demos/sp4/run.ipynb) ([Colab](https://colab.research.google.com/drive/1kYudJgbuhGBrhTTFs_HzT68LxFcVkJPu))
   * [Standard Problem #5](demos/sp5/run.ipynb) ([Colab](https://colab.research.google.com/drive/1RXlrHUtB39aHtyp2btk3GNEBS0f5ZDFk))
   * [Standard Problem Domainwall Pinning](demos/sp_domainwall_pinning/run.ipynb) ([Colab](https://colab.research.google.com/drive/1LgIX3o4e_6bww-RtIzJLX38QabUC5QMB))
   * [Standard Problem DMI](demos/sp_DMI/run.ipynb) ([Colab](https://colab.research.google.com/drive/1-5KuQ9GB3UeIfw4hCBN58fj2NvlXD28W))
   * [Standard Problem FMR](demos/sp_FMR/run.ipynb) ([Colab](https://colab.research.google.com/drive/1mN56sxjhgPuLA5yB7z3skmZ2cy733BbS))
   * [Standard Problem RKKY](demos/rkky/run.ipynb) ([Colab](https://colab.research.google.com/drive/1SIdiiz8plOI0SG3HhxNJYOxbknG178Qo))
   * [Stochastic Integration](demos/langevin/run.ipynb) ([Colab](https://colab.research.google.com/drive/1RlDaxgjqrZzerBFffDL7lQJHtEOm6v0q))
   * [Dispersion Calculator](demos/dispersion_calculator.ipynb) ([Colab](https://colab.research.google.com/drive/1B3sSPnm_Nycbka_Fa54INtXD2nZr8Mb2))

Installation
------------
For a clean and independent system, we start with a clean virtual python environment (this step could be omitted, if you would like to install magnum.np into the global python environment)

    mkdir venv
    python -m venv venv
    source venv/bin/activate


### from Python Package Index (PyPi)
In order to install a release versions of magnum.np one simple uses:

    pip install magnumnp

You can also easily install different versions from private repositories. E.g. use the following command to install the latest version of the main branch:

    pip install git+https://gitlab.com/magnum.np/magnum.np@main


### from source code (gitlab.com)
More advanced users can also install magnum.np from source code.
It can be downloaded from https://gitlab.com/magnum.np/magnum.np .

After activating the virtual environment magnum.np can be installed using the pip -e option which allows to easily modify the source code:

    pip install -e .

Note that a default version of [pytorch](http://www.pytorch.org) is included in magnum.np's dependecy list. If you would like to uses a specific pytorch version (fitting your installed CUDA library) it needs to be installed in advance.

### run remotely via Google Colab
---------------------------------
Magnum.np could also be used without any hardware by executing it remotely on resources provided by [Google Colab](https://drive.google.com/drive/folders/1Ymvx9bi0qQqW-zlOws0ahFJqoE3JCFd9?usp=share_link). The platform offers different runtime types like CPU(None), GPU or TPU. This allows users to directly test magnum.np, whithout needing their own hardware. Advanced users can use Google Colab(Pro), which provides access to current GPUs like the A100.

Some jupyter-notebook examples are included in the [demo](demos/README.md) directory, which also include links to Colab, where they can directly be run without any local installation.


Example
-------
The following demo code shows the solution of the muMAG Standard Problem #5 and can be found in the demos directory:

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
For example the following commands build an HTML documentation of the actual source code and stores it in the `public` folder:

    sphinx-build -b html docs public

Alternatively, the latest version of the documentation is always available on [https://magnum.np.gitlab.io/magnum.np/](https://magnum.np.gitlab.io/magnum.np/)


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
