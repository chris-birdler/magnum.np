.. image:: ./logo.png
  :width: 400
  :alt: magnum.np Logo

#####################################
magnum.np 1.1.1
#####################################

magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference method. It implements state-of-the-art algorithms and is based on `pytorch <http://www.pytorch.org/>`__, which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in Python which leads to very readable yet flexible code. Due to `pytorch <http://www.pytorch.org/>`__ integration, extensive postprocessing can be done directly in the simulations scripts. Alternatively, results can be written to PVD files and postprocessed with `Paraview <http://www.paraview.org/>`__. Furthermore `pytorch <http://www.paraview.org/>`__'s autograd feature makes it possible to solve inverse problems without significant modifications of the code. This manual is meant to give you both a quick start and a reference to magnum.np.

********
Features
********
* Explicit / Implicit time-integration of the Landau-Lifshitz-Gilbert Equation
* Fast FFT Demagnetization-field computation optimized for small memory footprint
* Fast FFT Oersted-field optimized for small memory footprint
* Periodic Boundary Conditions in 1D, 2D, and 3D (True and Pseudo-Periodic)
* Non-Equidistant Mesh for Multilayer Structures
* Arbitrary Material Parameters varying in space and time
* Spin-torque model by Zhang and Li
* Spin-Orbit torque (SOT)
* Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface, bulk, D2d)
* String method for energy barrier computations
* Sophisticated domain handling, e.g. for spatially varying material parameters
* Seamless VTK import / export via `pyvista <https://docs.pyvista.org/>`__
* Inverse Problems via `pytorch <www.pytorch.org/>`__'s autograd feature


*************
List of Demos
*************
* `RKKY <demos/rkky.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1HZyMxhG1HZCMsbdOohNqfmR0WUnoQ-ux>`__)
* `Softmagnetic Composite <demos/softmagnetic_composite.ipynb>`__ (`Colab <https://colab.research.google.com/drive/11dP3VrckM_hc24jP0sHIM_0_MlkNioRV#scrollTo=54df1504>`__)
* `Spin Orbit Torque <demos/sot.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1vfLhEUMGFCfJ-CB1QNKgZtLP6vkjf3hw#scrollTo=54df1504>`__)
* `Standard Problem #4 <demos/sp4.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1UNdTe5y41k_6HrZ7WsNASQSs5viUMzss>`__)
* `Standard Problem #5 <demos/sp5.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1AR7ksZUbThvJAn3mTgFKh6ka7dgXbUOb>`__)
* `Standard Problem Domainwall Pinning <demos/sp_domainwall_pinning.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1tMKEFtOfg7NSCGYUONeHLa0nC0AXVGRP#scrollTo=755d81e8>`__)
* `Standard Problem FMR <demos/sp_FMR.ipynb>`__ (`Colab <https://colab.research.google.com/drive/160QXXUkDEpd-GvZyI2PhLIdnSIDRVyG0>`__)


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   getting_started
   field_terms
   demos
   state_and_materials
   logging
   inverse_cube






********
Citation
********
If you use magnum.np in your work or publication, please cite the following reference:

[1] Bruckner, Florian, et al. "magnum.np -- A pytorch based GPU enhanced Finite Difference Micromagnetic Simulation Framework for High Level Development and Inverse Design", to be published (2023).


************
Contributing
************
Contributions are gratefully accepted.
The source code is hosted on `www.gitlab.com/magnum.np/magnum.np <http://www.gitlab.com/magnum.np/magnum.np>`__.
If you have any issues or question, just open an issue via gitlab.com.
To contribute code, fork our repository on gitlab.com and create a corresponding merge request.


******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
