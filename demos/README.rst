############
Introduction
############

The following demos should demonstrate the capabilities of magnum.np and also serve as a reference for the most relevant use cases. There are several ways to run the demos:


----------------------------
Run locally via command line
----------------------------
After downloading the source code and installing magnum.np in a virtual environment (or globally) the following commands can be used to run the demo from commandline

    python run.py
    python plot.py

The `run.py` script runs the actual simulation, whereas `plot.py` creates some output files which visualize the results of the simulation.


--------------------------------
Run locally via jupyter-notebook
--------------------------------
Jupyter notebooks allow to merge simulation and post-processing scripts in one file and add some additional documentation. The jupyter notebook can be openend locally via the following command:

    jupyter-notebook run.ipynb


-----------------------------
Run remotely via Google Colab
-----------------------------
The same notebook can be executed remotely on resources provided by Google Colab. The platform offers different runtime types like CPU(None), GPU or TPU. This allows users to directly test magnum.np, whithout needing their own hardware. Advanced users can use Google Colab(Pro), which provides access to current GPUs like the A100.


*************
List of Demos
*************

The following demos are also stored on Google Colab, where they can directly be run without any local installation:

* `Slonczewski Spin Torque 1 <https://florian98765.gitlab.io/magnum.np/notebooks/slonczewski1.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/slonczewski1.ipynb>`__)
* `Slonczewski Spin Torque 2 <https://florian98765.gitlab.io/magnum.np/notebooks/slonczewski2.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/slonczewski2.ipynb>`__)
* `Softmagnetic Composite <https://florian98765.gitlab.io/magnum.np/notebooks/softmagnetic_composite.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/softmagnetic_composite.ipynb>`__)
* `Spin Orbit Torque <https://florian98765.gitlab.io/magnum.np/notebooks/sot.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sot.ipynb>`__)
* `Standard Problem #4 <https://florian98765.gitlab.io/magnum.np/notebooks/sp4.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp4.ipynb>`__)
* `Standard Problem #5 <https://florian98765.gitlab.io/magnum.np/notebooks/sp5.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp5.ipynb>`__)
* `Standard Problem DMI <https://florian98765.gitlab.io/magnum.np/notebooks/sp_DMI.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_DMI.ipynb>`__)
* `Standard Problem Domainwall Pinning <https://florian98765.gitlab.io/magnum.np/notebooks/sp_domainwall_pinning.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_domainwall_pinning.ipynb>`__)
* `Standard Problem FMR <https://florian98765.gitlab.io/magnum.np/notebooks/sp_FMR.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/sp_FMR.ipynb>`__)
* `Standard Problem RKKY <https://florian98765.gitlab.io/magnum.np/notebooks/rkky.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky.ipynb>`__)
* `Stochastic Integration <https://florian98765.gitlab.io/magnum.np/notebooks/langevin.html>`__ (`Colab <https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/rkky.ipynb>`__)
* `Dispersion Calculator <demos/dispersion_calculator.ipynb>`__  (`Colab <https://colab.research.google.com/drive/1B3sSPnm_Nycbka_Fa54INtXD2nZr8Mb2>`__)
