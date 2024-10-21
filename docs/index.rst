.. image:: ./logo.png
  :width: 400
  :alt: magnum.np Logo

#####################################
magnum.np 2.0.1
#####################################

magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference
method. It implements state-of-the-art algorithms and is based on `pytorch <http://www.pytorch.org/>`__,
which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in
Python which leads to very readable yet flexible code. Due to `pytorch <http://www.pytorch.org/>`__
integration, extensive postprocessing can be done directly in the simulations scripts. Alternatively,
results can be written to PVD files and postprocessed with `Paraview <http://www.paraview.org/>`__.
Furthermore `pytorch <http://www.paraview.org/>`__'s autograd feature makes it possible to solve
inverse problems without significant modifications of the code. This manual is meant to give you
both a quick start and a reference to magnum.np.

**Version 2.0:** The magnum.np interface slightly changed since version 2.0.
Find more details about the necessary changes and the motivation `here <changes.rst>`__.
The following table summarizes the most important syntax changes:


********
Features
********
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
* Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface, bulk, D2d)
* String method for energy barrier computations
* `Eigenmode Solver <eigensolver.rst>`__ for efficient calculation of normal modes
* Sophisticated domain handling, e.g. for spatially varying material parameters
* efficient `Voronoi Code <voronoi.rst>`__ for 2D and 3D problems (including intergrain phase)
* Seamless VTK import / export via `pyvista <https://docs.pyvista.org/>`__
* Inverse Problems via `pytorch <www.pytorch.org/>`__'s autograd feature


*************
List of Demos
*************
Demo scripts for various applications are available in the [demo](demos/README.md) directory.

The following demos are also stored on Google Colab, where they can directly be run without any local installation:

* `Slonczewski Spin Torque <demos/slonczewski.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1KJCd-mnTaveZruLV37OWzYhtuDGuu3rH>`__)
* `Softmagnetic Composite <demos/softmagnetic_composite.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1HazB7ydSYZKbtrQoPc9xE3U0d7uc-1Ir>`__)
* `Spin Orbit Torque <demos/sot.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1OWMH0_qqxM73rB5gK5pi7nFRtO4nO_N8>`__)
* `Standard Problem #4 <demos/sp4.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1kYudJgbuhGBrhTTFs_HzT68LxFcVkJPu>`__)
* `Standard Problem #5 <demos/sp5.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1RXlrHUtB39aHtyp2btk3GNEBS0f5ZDFk>`__)
* `Standard Problem DMI <demos/sp_DMI.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1-5KuQ9GB3UeIfw4hCBN58fj2NvlXD28W>`__)
* `Standard Problem Domainwall Pinning <demos/sp_domainwall_pinning.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1LgIX3o4e_6bww-RtIzJLX38QabUC5QMB>`__)
* `Standard Problem FMR <demos/sp_FMR.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1mN56sxjhgPuLA5yB7z3skmZ2cy733BbS>`__)
* `Standard Problem RKKY <demos/rkky.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1SIdiiz8plOI0SG3HhxNJYOxbknG178Qo>`__)
* `Stochastic Integration <demos/langevin.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1RlDaxgjqrZzerBFffDL7lQJHtEOm6v0q>`__)
* `Dispersion Calculator <demos/dispersion_calculator.ipynb>`__  (`Colab <https://colab.research.google.com/drive/1B3sSPnm_Nycbka_Fa54INtXD2nZr8Mb2>`__)


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   changes
   getting_started
   field_terms
   demos
   state_and_materials
   eigensolver
   voronoi
   nonequidistant
   inverse_cube
   logging


********
Citation
********
If you use magnum.np in your work or publication, please cite the following reference:

[1] Bruckner, Florian, et al. "magnum.np -- A pytorch based GPU enhanced Finite Difference Micromagnetic Simulation Framework for High Level Development and Inverse Design", Scientific Reports volume 13, 12054 (2023).


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
