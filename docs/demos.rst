:tocdepth: 1

###############
Demos
###############

List of Demos
=============

The following demos should demonstrate the capabilities of magnum.np and also serve as a reference for the most relevant use cases:

Standard Problems
-----------------
* `Slonczewski Spin Torque 1           <https://github.com/magnumnp/magnumnp_demos/blob/main/slonczewski1.ipynb>`__           (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/slonczewski1.ipynb>`__)
* `Slonczewski Spin Torque 2           <https://github.com/magnumnp/magnumnp_demos/blob/main/slonczewski2.ipynb>`__           (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/slonczewski2.ipynb>`__)
* `Softmagnetic Composite              <https://github.com/magnumnp/magnumnp_demos/blob/main/softmagnetic_composite.ipynb>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/softmagnetic_composite.ipynb>`__)
* `Spin Orbit Torque                   <https://github.com/magnumnp/magnumnp_demos/blob/main/sot.ipynb>`__                    (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sot.ipynb>`__)
* `Standard Problem #4                 <https://github.com/magnumnp/magnumnp_demos/blob/main/sp4.ipynb>`__                    (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp4.ipynb>`__)
* `Standard Problem #5                 <https://github.com/magnumnp/magnumnp_demos/blob/main/sp5.ipynb>`__                    (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp5.ipynb>`__)
* `Standard Problem DMI                <https://github.com/magnumnp/magnumnp_demos/blob/main/sp_DMI.ipynb>`__                 (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_DMI.ipynb>`__)
* `Standard Problem Domainwall Pinning <https://github.com/magnumnp/magnumnp_demos/blob/main/sp_domainwall_pinning.ipynb>`__  (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_domainwall_pinning.ipynb>`__)
* `Standard Problem FMR                <https://github.com/magnumnp/magnumnp_demos/blob/main/sp_FMR.ipynb>`__                 (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_FMR.ipynb>`__)
* `Standard Problem RKKY               <https://github.com/magnumnp/magnumnp_demos/blob/main/rkky.ipynb>`__                   (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky.ipynb>`__)
* `Standard Problem biquadratic RKKY 1 <https://github.com/magnumnp/magnumnp_demos/blob/main/rkky_biquadratic1.ipynb>`__      (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky_biquadratic1.ipynb>`__)
* `Standard Problem biquadratic RKKY 2 <https://github.com/magnumnp/magnumnp_demos/blob/main/rkky_biquadratic2.ipynb>`__      (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky_biquadratic2.ipynb>`__)
* `Stochastic Integration              <https://github.com/magnumnp/magnumnp_demos/blob/main/langevin.ipynb>`__               (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky.ipynb>`__)
* `Dispersion Calculator               <https://github.com/magnumnp/magnumnp_demos/blob/main/dispersion_calculator.ipynb>`__  (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/dispersion_calculator.ipynb>`__)

Inverse Problems
----------------
* `Inverse Magnetization Reconstruction <https://github.com/magnumnp/magnumnp_demos/blob/main/inverse_magnetization_reconstruction.ipynb>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/inverse_magnetization_reconstruction.ipynb>`__)


Run the Demos
==============

There are several ways to run the demos:

Run locally via command line
----------------------------
After downloading the source code and installing magnum.np in a virtual environment (or globally) the following commands can be used to run the demo from commandline

.. code-block:: bash

    python run.py
    python plot.py

The `run.py` script runs the actual simulation, whereas `plot.py` creates some output files which visualize the results of the simulation.


Run locally via jupyter-notebook
--------------------------------
Jupyter notebooks allow to merge simulation and post-processing scripts in one file and add some additional documentation. The jupyter notebook can be openend locally via the following command:

.. code-block:: bash

    jupyter-notebook run.ipynb


Run remotely via Google Colab
-----------------------------
The same notebook can be executed remotely on resources provided by Google Colab. The platform offers different runtime types like CPU(None), GPU or TPU. This allows users to directly test magnum.np, whithout needing their own hardware. Advanced users can use Google Colab(Pro), which provides access to current GPUs like the A100.      
