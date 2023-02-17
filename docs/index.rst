#####################################
magnum.np 1.0.5
#####################################

magnum.np is a Python library for the solution of micromagnetic problems with the finite-difference method. It implements state-of-the-art algorithms and is based on `pytorch <http://www.pytorch.org/>`__, which allows to seamlessly run code either on GPU or on CPU. Simulation scripts are written in Python which leads to very readable yet flexible code. Due to `pytorch <http://www.pytorch.org/>`__ integration, extensive postprocessing can be done directly in the simulations scripts. Alternatively, results can be written to PVD files and postprocessed with `Paraview <http://www.paraview.org/>`__. Furthermore `pytorch <http://www.paraview.org/>`__'s autograd feature makes it possible to solve inverse problems without significant modifications of the code. This manual is meant to give you both a quick start and a reference to magnum.np.

********
Features
********
* Explicit / Implicit time-integration of the Landau-Lifshitz-Gilbert Equation
* Fast FFT Demagnetization-field computation optimized for small memory footprint
* Fast FFT Oersted-field optimized for small memory footprint
* Arbitrary Material Parameters variing in space and time
* Spin-torque model by Zhang and Li, Slonczewski
* Antiferromagnetic coupling layers (RKKY)
* Dzyaloshinskii-Moriya interaction (interface, bulk, D2d)
* String method for energy barrier computations
* Sophisticated domain handling, e.g. for spatially varying material parameters
* Seemingless VTK import / export via `pyvista <https://docs.pyvista.org/>`__
* Inverse Problems via `pytorch <www.pytorch.org/>`__'s autograd feature

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   getting_started
   field_terms


********
Citation
********
If you use magnum.np in your work or publication, please cite the following reference:

[1] Bruckner, Florian, et al. "magnum.np -- A pytorch based GPU enhanced Finite Difference Micromagnetic Simulation Framework for High Level Development and Inverse Design", to be published (2023).


************
Contributing
************
Contributions are gratefully accepted.
The source code is hosted on [www.gitlab.com/magnum.np/magnum.np](www.gitlab.com/magnum.np/magnum.np).
If you have any issues or question, just open an issue via gitlab.com.
To contribute code, fork our repository on gitlab.com and create a corresponding merge request.


******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
