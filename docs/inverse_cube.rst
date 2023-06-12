:tocdepth: 1

##################
Demo: Inverse Cube
##################

Usually the problems that are solved in micromagnetics are such that the final state is computed for given starting conditions, however, what if you wanted to find the necessary starting conditions for a certain state to be induced? The following demo code shows how magnum.np can be used to find the needed starting conditions for a given final state. 

The Code
********

First, all the necessary packages are imported. These are magnum.np and PyTorch which is a deep learning framework. The geometric function *sine* and *cosine* are imported separately from torch for readability. 

.. code-block:: python

  from magnumnp import *
  import torch
  from torch import sin, cos
  
A timer is set up in order to measure the execution time of the code. Set up a timer to measure the execution time of different parts of the code.

.. code-block:: python

  Timer.enable()
  
Next the mesh of the magnetic system is define. *n* defines the size of the magnetic material in terms of the number of cells in 3 dimensions and *dx* denotes the discretization. 

.. code-block:: python

  n  = (10, 10, 10)
  dx = (5e-9, 5e-9, 5e-9)
  mesh = Mesh(n, dx)
  
Simulation parameters are defined. These are the strength of the external magnetic field *h*, the relaxation time *tf* and the time step *dt*.  
  
.. code-block:: python  
  
  h = 300e-3/constants.mu_0
  tf = 1e-9
  dt = 1e-12
  
The initial state of the system is initialized. The *state.material* dictionary specifies parameters of the magnetic material.  
  
.. code-block:: python   

  state = State(mesh)
  state.material = {
      "Ms": 8e5,
      "A": 1.3e-11,
      "alpha": 1,
      "Ku": 1e5,
      "Ku_axis": state.Tensor([0,0,1])
      }
  
An initial magnetization *m0* is given. It is set to point in the positive z-direction, with a small perturbation in the x-direction.

.. code-block:: python

  m0 = state.Constant([0,0,0])
  m0[:, :, :, 2] = 1
  m0[:n[0]//2, :, :, 0] = 0.1
  m0[n[0]//2:, :, :, 0] = -0.1
  m0.normalize()
  
The desired final magnetization state *m_t* is set here as (1,0,1). This is the state for which optimal starting conditions are to be found.
  
.. code-block:: python

  m_t = state.Constant([1,0,1])
  m_t.normalize()
  
Both of these states are normalized.
  
Then, initial orientation angles *phi* and *theta* of the external field are given, *phi* being the polar angle and *theta* the azimuthal angle from spherical coordinates.
  
.. code-block:: python

  state.phi = state.Tensor([1.5708], requires_grad = True)
  state.theta = state.Tensor([1.5708], requires_grad = True)
  
The angles are defined as trainable variables by setting *requires_grad=True*. This enables automatic differentiation of these variables allowing them to be optimized in order to find the initial conditions for the desired final magnetization state. By using *requires_grad=True* all the operations that involve these tensors are recorded into a computational graph. Without this setting, PyTorch would not track the gradients of the loss with respect to these variables, which would not be updated during optimization. The loss quantifies how well the model's predictions fit the training data. The goal is of course to minimize the loss and thus find the optimal starting conditions. 

The external field is calculated with the lambda function *h_ext* using the angles *phi*, *theta* and the magnitude *h*. It takes the time *t* as input in case a time dependent external magnetic field is used. Here the field is constant over time. *h_ext* is returned in cartesian coordinates by converting the spherical coordinites.

.. code-block:: python  
  
  h_ext = lambda t: torch.stack([h*sin(state.theta)*cos(state.phi),
                                 h*sin(state.theta)*sin(state.phi),
                                 h*cos(state.theta)], dim=-1)
                                 
Then, the demagnetization, exchange, and anisotropy fields are defined. They describe the physical interactions between magnetic moments in the system. For more information about these fields please refer to :ref:`Sphinx Fields`. 

.. code-block:: python

  demag    = DemagField()
  exchange = ExchangeField()
  aniso    = UniaxialAnisotropyField()
  
The system is relaxed to the sate with the lowest energy using the LLGSolver.relax() method. This method solves the Landau-Lifshitz-Gilbert equation of motion for the magnetization of the system. During the relaxation process, m0 is updated to the final state, which is saved to a VTK file for visualization.

.. code-block:: python

  with torch.no_grad():
      state.m = m0
      llg = LLGSolver([demag, exchange, aniso])
      llg.relax(state)
      m0 = state.m
      write_vti(m0, "data/m0.vti")
  
After this follows the optimization process. This process is a loop in which the values of *phi* and *theta* are updated to minimize the loss between the final magnetization state that was just computed and the desired magnetization state *m_t*. 

First, an instance of the ExternalField class with the h_ext function as an argument is created. This allows the external field to be updated during the optimization loop to steer the system towards the desired final magnetization state. Furthermore, the LLGSolver object is created with a list of fields that include demagnetization, exchange, anisotropy, and external fields. The external field is defined by h_ext, which depends on the values of phi and theta. The solver object is created with the TorchDiffEqAdjoint solver, which is a solver that uses PyTorch's automatic differentiation capabilities to compute the gradients needed for optimization.

.. code-block:: python

  external = ExternalField(h_ext)
  llg = LLGSolver([demag, exchange, aniso, external], solver = TorchDiffEqAdjoint, adjoint_parameters = (state.phi, state.theta))

Three loggers: "phi", "theta", and "loss" are created to track the values of the orientation angles and the loss function during training, and two additional loggers "phi_grad" and "theta_grad" to track the gradients of the loss with respect to the orientation angles. The "phi_grad" and "theta_grad" loggers are defined using lambda functions that take in the current state of the system (which includes the orientation angles and their gradients) and return the gradient of the loss with respect to the corresponding angle. The ScalarLogger class is used to write the logged values. 

.. code-block:: python

  phi_grad = ("phi_grad", lambda state: state.phi.grad)
  theta_grad = ("theta_grad", lambda state: state.theta.grad)
  epoch_logger = ScalarLogger("data/epoch.dat", ["phi", "theta", "loss", phi_grad, theta_grad])

The optimization loop uses the Adam optimizer from PyTorch and the L1 loss function.. Therefore, an instance of *torch.optim.Adam* is created with *phi* and *theta* as the parameters to optimize. The learning rate is set to 0.05. *my_loss = torch.nn.L1Loss()* creates an instance of the loss function, which will be used to compute the loss between the final magnetization state *state.m* and the desired magnetization state *m_t* during the optimization loop. 

.. code-block:: python

  optimizer = torch.optim.Adam([state.phi, state.theta], lr=0.05)
  my_loss = torch.nn.L1Loss()
  
Next the omptimizatin loop is entered. The loop runs over 100 epochs. In each cycle the optimizer is zeroed and the initial magnetization state *m0* is set as the state of the system. The *state.t* variable is set to 0 to start the simulation from the beginning. 

.. cdoe-block:: python

  for epoch in range(100):
      print("epoch: ", epoch)
  
      optimizer.zero_grad()
      state.m = m0
      state.t = 0.

A logger is set up to record the time evolution of the magnetization during the optimization loop. The while loop updates the state of the system using the LLGSolver *llg.step(state, dt)* until the final time *tf* is reached. During each call, the magnetization state is evolved using the Landau-Lifshitz-Gilbert (LLG) equation and updated. *Timer("forward")* is simply used to measure the time it takes to run this loop.

.. code-block:: python

  logger = ScalarLogger("data/m_%03d.dat" % epoch, ['t', 'm'])
  with Timer("forward"):
      while state.t < tf:
          llg.step(state, dt)
          logger << state

*Timer("backward")* measures the time it takes to execute the following code. In thise code the loss between the final magnetization state of the system and the desired magnetization state *m_t* is calculated. The gradients of the loss with respect *phi* and *theta* are calculated by calling *state.loss.backward()*. Following this, the values of phi and theta are updated accordingly. Finally, the current epoch's data is saved to a log file using *epoch_logger << state*.

.. code-block:: python

  with Timer("backward"):
      state.loss = my_loss(state.m, m_t)
      state.loss.backward()

  optimizer.step()

  epoch_logger << state

Overall, the optimization process finds the values of *phi* and *theta* that minimize the difference between the final magnetization state of the system and the desired magnetization state, which would correspond to the optimal starting condition.


Lastly, the execution times of the different parts of the code are printed.

.. code-block:: python

  Timer.print_report()
  
The results of this simulation can be plotting using the following code:

.. code-block:: python

  import numpy as np
  import matplotlib.pyplot as plt
  
  data = np.loadtxt('data/epoch.dat')
  phi, theta, loss, phi_grad, theta_grad = data.T
  
  plt.plot(theta)
  plt.show()
  
Run the Simulation
******************

To run the simulation save the script to a file called *run.py* and enter the following into the command line

.. code-block:: python

  python run.py

See the Results
***************

After running run.py you can save the code for the plot in a file called plot.py. A plot of the results will be saved as *results.png*\ by entering the following in the command line:

.. code-block:: python

  python plot.py
