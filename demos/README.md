Introduction
============

The following demos should demonstrate the capabilities of magnum.np and also serve as a reference for the most relevant use cases. There are several ways to run the demos:


Run locally via command line
----------------------------
After downloading the source code and installing magnum.np in a virtual environment (or globally) the following commands can be used to run the demo from commandline

    python run.py
    python plot.py

The `run.py` script runs the actual simulation, whereas `plot.py` creates some output files which visualize the results of the simulation.


Run locally via jupyter-notebook
--------------------------------
Jupyter notebooks allow to merge simulation and post-processing scripts in one file and add some additional documentation. The jupyter notebook can be openend locally via the following command:

    jupyter-notebook run.ipynb


Run remotely via Google Colab
-----------------------------
The same notebook can be executed remotely on resources provided by Google Colab. The platform offers different runtime types like CPU(None), GPU or TPU. This allows users to directly test magnum.np, whithout needing their own hardware. Advanced users can use Google Colab(Pro), which provides access to current GPUs like the A100.



Documented Demos
================

The following demos are also stored on Google Colab, where they can directly be run without any local installation:

   * [RKKY](demos/rkky.ipynb) ([Colab] (https://colab.research.google.com/drive/1HZyMxhG1HZCMsbdOohNqfmR0WUnoQ-ux))
   * [Softmagnetic Composite](demos/softmagnetic_composite.ipynb) ([Colab](https://colab.research.google.com/drive/11dP3VrckM_hc24jP0sHIM_0_MlkNioRV))
   * [Spin Orbit Torque](demos/sot.ipynb) ([Colab](https://colab.research.google.com/drive/1vfLhEUMGFCfJ-CB1QNKgZtLP6vkjf3hw))
   * [Standard Problem #4](demos/sp4.ipynb) ([Colab](https://colab.research.google.com/drive/1UNdTe5y41k_6HrZ7WsNASQSs5viUMzss))
   * [Standard Problem #5](demos/sp5.ipynb) ([Colab](https://colab.research.google.com/drive/1AR7ksZUbThvJAn3mTgFKh6ka7dgXbUOb))
   * [Standard Problem Domainwall Pinning](demos/sp_domainwall_pinning.ipynb) ([Colab](https://colab.research.google.com/drive/1tMKEFtOfg7NSCGYUONeHLa0nC0AXVGRP))
   * [Standard Problem FMR](demos/sp_FMR.ipynb) ([Colab](https://colab.research.google.com/drive/160QXXUkDEpd-GvZyI2PhLIdnSIDRVyG0))
