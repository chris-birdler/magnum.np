Introduction
============

The following demos should demonstrate the capabilities of magnum.np and also serve as a reference for the most relevant use cases. There are several ways to run the demos:


Run locally via command line
----------------------------
After downloading the source code and installing magnum.np in a virtual environment (or globally) the following commands can be used to run the demo from commandline

    python run.py
    python plot.py

The `run.py` script runs the actual simulation, whereas `plot.py` creates some output files which visualize the results of the simulation.


Run locally via jupiter-notebook
--------------------------------
Jupiter notebooks allow to merge simulation and post-processing scripts in one file and add some additional documentation. The jupiter notebook can be openend locally via the following command:

    jupyter-notebook run.ipynb


Run remotely via Google Colab
-----------------------------
The same notebook can be executed remotely on resources provided by Google Colab. The platform offers different runtime types like CPU(None), GPU or TPU. This allows users to directly test magnum.np, whithout needing their own hardware. Advanced users can use Google Colab(Pro), which provides access to current GPUs like the A100.      



Documented Demos
================

[Standard Problem #4](sp4/run.ipynb)
