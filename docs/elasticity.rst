.. module:: magnumnp

:tocdepth: 1

#################
Linear Elasticity
#################

In ferromagnetic materials, the magnetic and elastic degrees of freedom are coupled by the (inverse) magnetostrictive effect. Due to this coupling, acoustic waves can excite spin waves in ferromagnetic materials. The following demo shows how the self-consistent dynamics of the magnetization and mechanical displacement can be obtained for an acoustic bulk mode.

******************
State and Material
******************

To investigate the coupling with a bulk mode, a quasi-1D mesh with periodic boundary conditions is sufficient.

.. code-block:: python

    mesh = Mesh(n, dx, pbc=(0,1,1))
    state = State(mesh)

The coupled magnetoelastic solver requires additional material parameters. These are the mass density ``"rho"`` given in :math:`\text{kg/m}^3`, the stiffness matrix ``"C"`` given in :math:`\text{N/m}^2` and written in Voigt notation, the dimensionless coupling parameters ``"lambda_100"`` and ``"lambda_111"`` as well as a phenomenological damping parameter ``"eta"``, which is set to zero in this example and is mainly used to relax structures.

The dynamics of the mechanical degrees of freedom further depends on the magnetization, and is also explicitly dependent on the exchange constant ``"A"`` at boundary and material interface cells.

.. code-block:: python
    state.material = {
            "alpha": state.Constant(alpha),
            "Ms": state.Constant(Ms),
            "A": state.Constant(Aex),
            "rho": state.Constant(rho),
            "C": state.Constant(C_cubic(C11, C12, C44)),
            "lambda_100": state.Constant(lambda_100),
            "lambda_111": state.Constant(lambda_111),
            "eta": state.Constant(eta)
        }

The mechanical degrees of freemdoms need to be initialized together with the magnetization as members of ``state``. The displacement is ``state.ud``, the momentum density is ``state.pd``. For this example, we do not relax the structure but instead compute the initial displacement due to the magnetic strain.

.. code-block:: python

    # initial conditions
    state.m = state.Constant((0.,0.,1.))
    state.pd = state.Constant((0.,0.,0.))

    # ... initial displacement due to magnetic strain
    kx_x = -0.5*lambda_100
    kx_y = -0.5*lambda_100
    kx_z = lambda_100
    kx = kx_x + (C12/C11)*(kx_y + kx_z)

    x = mesh.SpatialCoordinate()[0]
    state.ud[...,0] = -((Lx-0.5*dx[0])-x)*kx

*************************************
Field Terms and Solver Initialization
*************************************

The coupled solver automatically considers the magnetic state when updating the displacement and momentum density, but the magnetelastic field term needs to be added manually to the list of considered field terms by the user.

The parameter ``C_sym`` in the constructor of the coupled solver allows the user to specify which terms of the stiffness matrix are non-zero. Specifying this can result in better timings. We also set ``iteration_depth=0``. This indicates that gradient components of the displacement do not need to be computed iteratively, as it would be necessary for a 2D or 3D problem. The reason for this is that in our 1D example, only spatial derivatives of :math:`x` are non-zero, and no transversal gradient components can enter the jump conditions of the spatial derivatives.

.. code-block:: python

    bias = ExternalField(state.Constant((0, 0, h0))
    exchange = ExchangeField()
    magEl = MagnetoElasticField()
    h_terms = [bias, exchange, magEl]

    llg = LLGWithLESolver(h_terms, C_sym="cubic", iteration_depth=0, boundary_nodes=3, dt=dt)

*******************************************
Boundary Conditions for the Elastic Problem
*******************************************

When we initialized the coupled solver, we specified the number of boundary nodes used for boundary value computations by setting ``boundary_nodes = 3``. The number of boundary nodes used (counting inwards from the boundary) can be set between 1 and 3. Using three nodes will result in the highest accuracy, but it also requires that the first three cells from all boundaries belong to the same material.

| When ``boundary_nodes = 3`` is set, the boundary node values of the strain i.e. the spatial derivatives of the displacement components will be calculated using a second order scheme. For open surfaces, the boundary node values of the force density will be updated using these boundary node values, the strain on the neighbouring bulk node and the traction set as Neumann boundary conditions on the surface.
| When ``boundary_nodes = 2`` is set, the strain is calculated from forward and backward differences formulas at the boundary nodes. At open surfaces, the force density is calculated from the traction given by the Neumann boundary condition and the strain on the first bulk node. Therefore, the boundary values of the strain, that are computed from forward or backward differences, are skipped when computing the force density on the boundary nodes.
| The default is ``boundary_nodes = 1``, where forward and backward differences are used for the boundary values of the displacement gradients and the force density. This method is inaccurate, but applicable even if the first two layers of cells counting from the surface belong to different materials.

This setting informs the solver how to handle the boundary conditions, but the boundary conditions themselves need to be added in form of dictionary as member ``state.bcs``. The keys specify the type of the boundary condition. The key ``"ud"`` is associated with Dirichlet boundary conditions, and the corresponding value has to be a list of ``BC``. The key ``"t"`` is associated to Neumann boundary conditions and the corresponding value needs to be a list of ``PlaneBC``. The main difference between the two is *where* the boundary conditions are defined. ``BC`` is used to set value of ``state.ud`` on the cell center while ``PlaneBC`` is used to set the traction at a cell interface.

.. csv-table:: Boundary Conditions
   :header: "Type", "key", "Type", defined at

   "Dirichlet BCs", ``"ud"``, ":code:`list[BC]`", "cell center"
   "Neumann BCs", ``"t"``, ":code:`list[PlaneBC]`", "cell interface"


.. code-block:: python

    u_mask = torch.zeros(n, dtype=torch.bool)
    u_mask[-1] = 1
    u_cond = lambda state : torch.zeros((n[1], n[2], 3))

    dirichlet_bcs = []
    dirichlet_bcs.append(BC(u_mask, u_cond))

    # ... left shear traction
    def t0_cond(state):
        nx, ny, nz = state.mesh.n

        t0 = u0 * C44 * k

        t_in = torch.zeros((ny, nz, 3))
        t_in[:,:,2] = t0 * np.sin(float(state.t) * omega)

        return t_in

    neumann_bcs = []
    neumann_bcs.append(PlaneBC(Plane(0, 0, -1), t0_cond))

    state.bcs = dict()
    state.bcs["ud"] = dirichlet_bcs
    state.bcs["t"] = neumann_bcs

.. note::
    If no Neumann boundary condition and no Dirichlet boundary condition is set at a boundary cell, homogeneous Neumann boundary conditions are set automatically.

.. note::
    It is also possible to introduce an airbox, or individual cells of air (cells or layers in which the stiffness matrix is set to zero) to enforce homogeneous Neumann boundary conditions. This method generally yields the most stable results, since the boundary conditions are enforced already from the jump conditions when calculating the spatial derivatives.


.. autoclass:: BC

.. autoclass:: PlaneBC

.. autoclass:: Plane


*********************
Time Step and Logging
*********************

Time steps are performed analogous to time stepping with ``LLGSolver``.

.. code-block:: python

    llg.step(state, dt)

Besides the relative tolerance ``rtol``, the absolute tolerance can be specified for the magnetization, displacement and momentum density individually by the arguments ``atol_m``, ``atol_ud``, and ``atol_pd``, respectively.

For logging purposes, the coupled solver allows to compute multiple energy terms from builtin methods.

.. csv-table:: Energy terms
   :header: "Method", "Energy Terms"

   ":code:`LLGWithLESolver.E(state)`", "Sum of the energy of the magnetic field terms"
   ":code:`LLGWithLESolver.U_el(state)`", "Elastic potential energy"
   ":code:`LLGWithLESolver.U(state)`", "Mechanical potential energy"
   ":code:`LLGWithLESolver.T_el(state)`", "Kinetic Energy"

*************
Complete Code
*************

The complete code can be viewed here: :download:`run.py <../demos/linear_elasticity/bulk_shear_mode/run.py>`.

##############
Coupled Solver
##############

.. autoclass:: LLGWithLESolver
