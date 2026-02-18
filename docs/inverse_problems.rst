.. currentmodule:: magnumnp

################
Inverse Problems
################

In classical (forward) micromagnetic simulations the effective fields and the resulting magnetization dynamics are computed for given initial conditions and material parameters.
*Inverse problems* reverse this process: given some target observable — a desired final state, a measured stray field, or an experimentally obtained spectrum — the goal is to find the unknown inputs (initial conditions, material constants, applied-field protocols, etc.) that reproduce it.

Because every field computation in *magnum.np* is implemented as a differentiable PyTorch operation, inverse problems can be solved by gradient-based optimization with minimal changes to a regular simulation script.
The following sections outline the general recipe; concrete worked-out examples are linked below.

Demos
-----

.. toctree::
   :maxdepth: 1

   Magnetization Reconstruction from Stray-Field Data <inverse_mag_rec>
   Optimal External-Field Angles <inverse_cube>
  


General Structure of an Inverse Problem
========================================

Every inverse problem solved with *magnum.np* follows the same five-step pattern.


1. Forward Problem
------------------

Set up the micromagnetic system (mesh, material, fields) and define the forward computation exactly as for a regular simulation.
This can be a time integration of the LLG equation, the evaluation of a demagnetization field, an energy calculation, or any other operation provided by *magnum.np*.

.. code-block:: python

   from magnumnp import *
   import torch

   mesh = Mesh(n, dx)
   state = State(mesh)
   state.material = { ... }

   demag    = DemagField()
   exchange = ExchangeField()
   # ... additional field terms

   def forward_model(params):
       # Run the micromagnetic simulation to compute the quantity of interest
       ...
       return simulated_quantity   # will be compared to a target later


2. Trainable Parameters (``requires_grad``)
--------------------------------------------

Identify which quantities are unknown and mark them as trainable by setting ``requires_grad=True``.
PyTorch will then record every operation involving these tensors into a computational graph, enabling automatic differentiation.

.. code-block:: python

   # Example: optimize external-field angles
   state.phi   = state.Tensor([1.5708], requires_grad=True)
   state.theta = state.Tensor([1.5708], requires_grad=True)

   # Example: optimize a full magnetization field
   m_opt = torch.zeros(nx, ny, nz, 3)
   m_opt.requires_grad = True

Any tensor can be made trainable — scalar parameters, spatially varying material constants, or the magnetization itself.


3. Loss Function
-----------------

Define a scalar *loss* that quantifies how far the current forward-model output is from the target.
Typical choices include the L1 or L2 distance between simulated and measured quantities, optionally combined with physics-based regularization terms (exchange energy, DMI energy, etc.).

.. code-block:: python

   criterion = torch.nn.L1Loss()
   loss = criterion(simulated_quantity, target_quantity)

   # Optional: add regularization
   loss = loss + lambda_reg * regularization


4. Backward Pass (``backward``)
---------------------------------

Calling ``loss.backward()`` triggers PyTorch's automatic differentiation (autograd) engine.
It traverses the recorded computational graph in reverse and computes the gradient of the loss with respect to every trainable parameter.
Because all *magnum.np* field computations are built from differentiable PyTorch primitives (FFTs, convolutions, stencil operations), gradients flow seamlessly through the entire micromagnetic forward model.

.. code-block:: python

   loss.backward()


5. Optimizer Step (``zero_grad`` / ``step``)
----------------------------------------------

A PyTorch optimizer updates the trainable parameters in the direction that reduces the loss.
Before each forward evaluation the optimizer's gradient buffers are cleared with ``zero_grad()``, and after the backward pass the parameters are updated with ``step()``.

.. code-block:: python

   optimizer = torch.optim.Adam([state.phi, state.theta], lr=0.05)

   for epoch in range(n_epochs):
       optimizer.zero_grad()

       # --- forward problem ---
       ...

       # --- loss ---
       loss = ...

       # --- backward ---
       loss.backward()

       # --- update ---
       optimizer.step()


Putting It All Together
========================

The complete optimization loop therefore reads:

.. code-block:: python

   from magnumnp import *
   import torch

   # 1. Forward problem setup
   mesh  = Mesh(n, dx)
   state = State(mesh)
   state.material = { ... }
   fields = [DemagField(), ExchangeField(), ...]

   # 2. Trainable parameters
   params = state.Tensor([...], requires_grad=True)

   # 3. Loss function
   criterion = torch.nn.L1Loss()

   # 5. Optimizer
   optimizer = torch.optim.Adam([params], lr=0.05)

   for epoch in range(n_epochs):
       optimizer.zero_grad()          # clear gradients

       # --- run forward model ---
       ...

       loss = criterion(output, target) # compute loss
       loss.backward()                  # compute gradients
       optimizer.step()                 # update parameters

Because *magnum.np* delegates all numerical work to PyTorch, any optimizer from ``torch.optim`` can be used (Adam, L-BFGS, SGD, ...) and advanced techniques such as learning-rate scheduling are available out of the box.
