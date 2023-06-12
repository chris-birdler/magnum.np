:tocdepth: 1

###############
Demos
###############

List of Demos
=============

The following demos should demonstrate the capabilities of magnum.np and also serve as a reference for the most relevant use cases:


* `RKKY <demos/rkky.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1SIdiiz8plOI0SG3HhxNJYOxbknG178Qo>`__)

 
* `Softmagnetic Composite <demos/softmagnetic_composite.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1HazB7ydSYZKbtrQoPc9xE3U0d7uc-1Ir>`__)
 
 
* `Spin Orbit Torque <demos/sot.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1OWMH0_qqxM73rB5gK5pi7nFRtO4nO_N8>`__)

 
* `Standard Problem #4 <demos/sp4.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1kYudJgbuhGBrhTTFs_HzT68LxFcVkJPu>`__)
 

 
* `Standard Problem #5 <demos/sp5.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1RXlrHUtB39aHtyp2btk3GNEBS0f5ZDFk>`__)

* `Standard Problem DMI <demos/sp_DMI.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1-5KuQ9GB3UeIfw4hCBN58fj2NvlXD28W>`__)
 
 
* `Standard Problem Domainwall Pinning <demos/sp_domainwall_pinning.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1LgIX3o4e_6bww-RtIzJLX38QabUC5QMB>`__)
 

* `Standard Problem FMR <demos/sp_FMR.ipynb>`__ (`Colab <https://colab.research.google.com/drive/1mN56sxjhgPuLA5yB7z3skmZ2cy733BbS>`__)



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


