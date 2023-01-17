# magnum.np

[![Documentation Status](https://readthedocs.com/projects/magnumnp-magnumnp/badge/?version=latest&token=14bfc008bc74536d868c875c61777aea0733cb782c08690a29aecf5761e21f1f)](https://magnumnp-magnumnp.readthedocs-hosted.com/en/latest/?badge=latest)
      
magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference method. It implements state-of-the-art algorithms and is based on [pytorch](http://www.pytorch.org/), which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in Python which leads to very readable yet flexible code. Due to [pytorch](http://www.pytorch.org/) integration, extensive postprocessing can be done directly in the simulations scripts. Alternatively, results can be written to PVD files and postprocessed with [Paraview](http://www.paraview.org/). Furthermore [pytorch](http://www.paraview.org/)'s autograd feature makes it possible to solve inverse problems without significant modifications of the code. This manual is meant to give you both a quick start and a reference to magnum.np.

********
Features
********
* Explicit time-integration of the Landau-Lifshitz-Gilbert Equation
* Implicit time-integration of the Landau-Lifshitz-Gilbert Equation (currently only on CPU, using scipy solvers)
* Fast FFT Demagnetization-field computation optimized for small memory footprint
* Fast FFT Oersted-field  optimized for small memory footprint
* Arbitrary field terms varying in space and time
* Spin-torque model by Zhang and Li
* Spin-torque model by Slonczewski
* Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface and bulk - Experimental!)
* String method for energy barrier computations
* Sophisticated domain handling, e.g. for spatially varying material parameters
