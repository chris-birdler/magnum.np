:tocdepth: 1

#######################################
Inverse Magnetization Reconstruction
#######################################

Magnetic inverse problems aim to infer unknown magnetization configurations from indirect measurements such as stray-field data.
While classical (forward) micromagnetic simulations compute the field produced by a known magnetization, the inverse problem asks:
*given stray-field measurements, what is the underlying magnetization?*

Thanks to PyTorch's automatic differentiation (autograd), *magnum.np* can solve such inverse problems by formulating them as optimization tasks.
The key idea is to treat the unknown magnetization as a set of trainable parameters and minimize a loss function that quantifies the mismatch between simulated and measured stray fields, optionally regularized by micromagnetic energy terms.

This page demonstrates the reconstruction of a 2D magnetization texture from a single stray-field projection, using a Cartesian parameterization with physics-based regularization.
For a simpler introductory example of inverse problems with *magnum.np* (finding optimal external-field angles), see :doc:`inverse_cube`.


Problem Statement
==================

Consider a thin magnetic film whose magnetization :math:`\mathbf{m}(\mathbf{r})` is unknown.
A sensor (e.g.\ a nitrogen-vacancy magnetometer) measures the projection of the stray field onto a fixed axis :math:`\mathbf{n}_\mathrm{NV}` at a measurement plane above the sample:

.. math::

   h_\mathrm{meas}(\mathbf{r}_\parallel) = \mathbf{H}^\mathrm{dem}(\mathbf{r}_\parallel, z_\mathrm{meas}) \cdot \mathbf{n}_\mathrm{NV}.

The reconstruction task is to find the magnetization :math:`\mathbf{m}` that minimizes

.. math::

   \mathcal{L}(\mathbf{m}) \;=\; \underbrace{\bigl\| h_\mathrm{sim}(\mathbf{m}) - h_\mathrm{meas} \bigr\|_1}_{\text{data fidelity}}
   \;+\; \lambda \underbrace{\bigl( E_\mathrm{exc} + E_\mathrm{DMI} + E_\mathrm{aniso} + E_\mathrm{dem} \bigr)}_{\text{physics regularization}},

subject to the unit-length constraint :math:`|\mathbf{m}| = 1` everywhere.

The data-fidelity term drives the solution towards consistency with the measurement, while the physics regularization steers the reconstruction towards micromagnetically plausible states and acts as a prior to resolve the inherent non-uniqueness of the inverse problem.


The Code
*********

Setup
-----

First, the necessary packages are imported: *magnum.np* for micromagnetic computations, PyTorch for automatic differentiation, and *pathlib* for path handling.

.. code-block:: python

   from magnumnp import *
   import torch
   import pathlib

   try:
       this_dir = pathlib.Path(__file__).resolve().parent
   except:
       this_dir = pathlib.Path().resolve()


Normalization Helper
---------------------

A helper function normalizes unconstrained 3-component vectors to the unit sphere. This is used throughout the optimization to enforce the constraint :math:`|\mathbf{m}| = 1`.

.. code-block:: python

   def normalize_to_unit(m):
       """
       Given unconstrained m of shape (Nx, Ny, Nz, 3), return unit-normalized magnetization.
       """
       return m / m.norm(dim=-1, keepdim=True)


Simulation Parameters
----------------------

The computational grid, material constants, and optimization hyperparameters are defined. The material parameters correspond to a thin film with perpendicular magnetic anisotropy and interfacial Dzyaloshinskii–Moriya interaction (DMI), which supports chiral domain structures such as Néel skyrmions.

.. code-block:: python

   # Grid
   dx, dy, dz = 3e-9, 3e-9, 2.4e-9

   # Material
   Ms = 1.4e6
   A = 1.5e-11
   Di = 3e-3
   Ku = 1.1e6
   Ku_axis = (0, 0, 1)

   # Optimization
   n_epochs = 30
   lr = 1.0
   max_iter = 100
   history_size = 100

   # Regularization weight (balances data fidelity vs. physics energy)
   lambda_reg = 5e17

The regularization weight :math:`\lambda` controls the trade-off between fitting the measurement data and satisfying micromagnetic energy constraints. A larger :math:`\lambda` enforces more physically plausible solutions at the cost of slightly worse data fit.


Loading the Ground Truth
--------------------------

A pre-relaxed micromagnetic state is loaded from a VTI file. In a real experiment this would be replaced by actual measurement data; here it serves as a synthetic ground truth to validate the reconstruction.

.. code-block:: python

   mesh, fields = read_vti(this_dir / "data" / "m_true.vti")
   m_true = fields["f000"]


Domain and Material Setup
---------------------------

The simulation domain is constructed as a 3D grid with a single material layer at the bottom and vacuum above it. The vacuum region is needed so that the demagnetization field can be evaluated at a measurement plane above the sample. Material parameters are set to near-zero in the vacuum and to the physical values inside the material domain.

.. code-block:: python

   nx, ny, nz_material = m_true.shape[:3]
   nz_vacuum = 30
   nz = nz_material + nz_vacuum

   mesh = Mesh((nx, ny, nz), (dx, dy, dz))
   state = State(mesh)

   # Active domain: bottom layer only
   domain = state.Constant(False, dtype=torch.bool)
   domain[:, :, 0] = True

   state.material = {
       "Ms": 1e-20,
       "A": 1e-20,
       "Di": 1e-20,
       "Ku": 1e-20,
       "Ku_axis": torch.tensor(Ku_axis).to(state.device),
   }

   state.material["A"][domain] = torch.tensor(A)
   state.material["Ms"][domain] = torch.tensor(Ms)
   state.material["Di"][domain] = torch.tensor(Di)
   state.material["Ku"][domain] = torch.tensor(Ku)

Note that the material parameters inside the vacuum are set to very small positive values (``1e-20``) rather than exactly zero, to avoid division-by-zero issues in the field computations.


Computing the Target Stray Field
----------------------------------

The ground-truth magnetization is placed into the material layer and the resulting demagnetization field is computed at the top of the simulation domain (the measurement plane). The field is projected onto the NV-center axis :math:`\mathbf{n}_\mathrm{NV}` to produce a scalar measurement map, and ``.detach()`` is called to exclude the target from the computational graph.

.. code-block:: python

   state.m = state.Constant([0.0, 0.0, 0.0])
   state.m[:, :, :1, :] = m_true.clone()

   # Normal vector for projection (NV-center axis)
   nv = torch.tensor([0.819152, 0, 0.573576])

   # Compute demag field at measurement plane
   demag_field = DemagField()
   h_demag = demag_field.h(state)
   h_vec = h_demag[:, :, -1, :]

   # Project onto NV axis
   h_measurement_target = (h_vec @ nv).unsqueeze(-1).detach()

The field objects used for physics-based regularization are also initialized here:

.. code-block:: python

   exchange_field = ExchangeField()
   dmi_field = InterfaceDMIField()
   aniso_field = UniaxialAnisotropyField()


Setting Up the Optimizer
--------------------------

The reconstruction starts from a uniform :math:`+z` initial guess. The magnetization is stored as an unconstrained tensor ``m_opt`` with ``requires_grad=True``, enabling PyTorch to track all operations for gradient computation.

The L-BFGS optimizer is used because it converges quickly for smooth, high-dimensional problems. It approximates the inverse Hessian using a history of past gradient evaluations, providing quasi-Newton convergence without the memory cost of storing the full Hessian matrix. The Strong Wolfe line search ensures sufficient descent at each step.

.. code-block:: python

   m_opt = torch.zeros(nx, ny, 1, 3)
   m_opt[..., 2] = 1.0

   with torch.no_grad():
       m_0 = normalize_to_unit(m_opt).clone()

   m_opt.requires_grad = True

   optimizer = torch.optim.LBFGS(
       [m_opt],
       lr=lr,
       max_iter=max_iter,
       history_size=history_size,
       line_search_fn="strong_wolfe",
   )


Closure Function
------------------

L-BFGS requires a *closure* — a callable that computes the loss and its gradients. Inside the closure:

1. The unconstrained ``m_opt`` is normalized to the unit sphere.
2. The normalized magnetization is placed into the simulation state.
3. The **data-fidelity loss** is computed as the L1 distance between the simulated and target stray-field projections.
4. The **physics regularization** sums the exchange, DMI, anisotropy, and demagnetization energies.
5. ``loss.backward()`` triggers PyTorch's autograd to compute :math:`\partial\mathcal{L}/\partial\tilde{\mathbf{m}}` through the entire micromagnetic forward model.

.. code-block:: python

   def closure():
       optimizer.zero_grad()

       m_norm = normalize_to_unit(m_opt)

       m_full = state.Constant([0.0, 0.0, 0.0])
       m_full[:, :, :1, :] = m_norm
       state.m = m_full

       # Data fidelity loss
       h_demag_current = demag_field.h(state)
       h_vec_current = h_demag_current[:, :, -1, :]
       h_proj_current = (h_vec_current @ nv).unsqueeze(-1)
       data_loss = torch.nn.functional.l1_loss(h_proj_current, h_measurement_target)

       # Physics energy regularization
       energy = (exchange_field.E(state)
                 + dmi_field.E(state)
                 + aniso_field.E(state)
                 + demag_field.E(state))

       loss = data_loss + lambda_reg * energy

       loss.backward()

       closure.data_loss = data_loss.item()
       closure.energy = energy.item()
       closure.total_loss = loss.item()
       return loss

Because all field computations in *magnum.np* are implemented as differentiable PyTorch operations, the gradients flow seamlessly through the FFT-based demagnetization, the exchange stencil, and all other field terms.


Running the Optimization
--------------------------

The optimization loop calls ``optimizer.step(closure)`` for each epoch. The L-BFGS optimizer internally evaluates the closure multiple times per step (controlled by ``max_iter``) to find the optimal step size along the descent direction.

.. code-block:: python

   print(f"{'Epoch':>6} | {'Total Loss':>12} | {'Data Loss':>12} | {'Energy':>12} | {'m_diff':>12}")
   for epoch in range(1, n_epochs + 1):
       optimizer.step(closure)

       # Snapshot after first optimizer step
       if epoch == 1:
           with torch.no_grad():
               m_rec_1 = normalize_to_unit(m_opt).clone()

       # Compute current reconstruction difference
       with torch.no_grad():
           m_rec_current = normalize_to_unit(m_opt)
           m_diff_current = (m_rec_current[...,0,:] - m_true[...,0,:]).norm(dim=-1).mean().item()

       if epoch == 1 or epoch % 10 == 0:
           print(f"{epoch:6d} | {closure.total_loss:12.4e} | {closure.data_loss:12.4e} | {closure.energy:12.4e} | {m_diff_current:12.4e}")


Evaluating and Saving Results
-------------------------------

After the optimization the final reconstruction is compared against the ground truth, and some intermediate states are saved for visualization.

.. code-block:: python

   with torch.no_grad():
       m_rec = normalize_to_unit(m_opt)

       m_diff = (m_rec - m_true).norm(dim=-1).mean().item()
       print(f"\nFinal Mean Vector Diff vs Truth: {m_diff:.4f}")

       h_final = (demag_field.h(state)[:, :, -1, :] @ nv).unsqueeze(-1)
       h_err = torch.abs(h_final - h_measurement_target).mean().item()
       print(f"Final H-field MAE: {h_err:.4e}")

       torch.save({
           "m_true": m_true,
           "m_0": m_0,
           "m_rec_1": m_rec_1,
           "m_rec": m_rec,
           "h_true": h_true,
           "h_0": h_0,
           "h_rec_1": h_rec_1,
           "h_rec": h_rec,
       }, this_dir / "data" / "results.pt")
       print("Saved results to data/")


Plotting the Results
---------------------

The reconstruction can be visualized with the following plotting script, which shows the three magnetization components and the stray-field projection for the ground truth, initial guess, first optimization step, and final reconstruction:

.. code-block:: python

   import torch
   import matplotlib.pyplot as plt

   data = torch.load("data/results.pt")
   m_true = data["m_true"].cpu()
   m_0 = data["m_0"].cpu()
   m_rec_1 = data["m_rec_1"].cpu()
   m_rec = data["m_rec"].cpu()
   h_true = data["h_true"].cpu()
   h_0 = data["h_0"].cpu()
   h_rec_1 = data["h_rec_1"].cpu()
   h_rec = data["h_rec"].cpu()

   comp_labels = ["$m_x$", "$m_y$", "$m_z$",
                   r"$\mathbf{H}^\mathrm{dem} \cdot \mathbf{n}_\mathrm{NV}$"]
   row_labels = ["Ground Truth", "Starting Guess", "After 1 Step", "Final"]
   row_data = [m_true, m_0, m_rec_1, m_rec]
   row_fields = [h_true, h_0, h_rec_1, h_rec]

   fig, axes = plt.subplots(4, 4, figsize=(20, 16))
   for row, (label, m, h) in enumerate(zip(row_labels, row_data, row_fields)):
       for i in range(3):
           im = axes[row, i].imshow(m[:, :, 0, i].T, origin="lower",
                                     cmap="RdBu_r", vmin=-1, vmax=1)
           if row == 0:
               axes[row, i].set_title(comp_labels[i], fontsize=12)
           if i == 0:
               axes[row, i].set_ylabel(label, fontsize=11, fontweight='bold')
           plt.colorbar(im, ax=axes[row, i], fraction=0.046)

       im = axes[row, 3].imshow(h.T, origin="lower", cmap="viridis")
       if row == 0:
           axes[row, 3].set_title(comp_labels[3], fontsize=12)
       plt.colorbar(im, ax=axes[row, 3], fraction=0.046)

   fig.suptitle("Magnetization Reconstruction Results", fontsize=16, fontweight='bold')
   fig.tight_layout()
   fig.savefig("data/results.png", dpi=150)


Run the Simulation
*******************

To run the reconstruction, save the script to a file called *run.py* and execute:

.. code-block:: bash

   python run.py


See the Results
****************

After running *run.py*, the plotting script can be executed to produce visualizations:

.. code-block:: bash

   python plot.py

This generates ``data/results.png`` showing the magnetization components and stray-field projections at each stage of the reconstruction, as well as ``data/error_map.png`` showing the pointwise reconstruction error.
