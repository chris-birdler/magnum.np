.. currentmodule:: magnumnp

:tocdepth: 1

#################
Linear Elasticity
#################

Magnetoelastic couplings are available via the :class:`.MagnetoElasticField`
and the :class:`.LinearMagnetoElasticField`.  The full field uses the strain
induced by the magnetization as well as the mechanical response, while the
linearized version only keeps the terms that are quadratic in the
magnetization.  Both fields require elastic material data (density and
stiffness tensor) and the magnetostriction constants
:math:`\lambda_{100}`/:math:`\lambda_{111}` to be present in
``state.material``.

Setting Up the Magnetoelastic State
===================================

The coupled solver expects the elastic degrees of freedom and material
parameters to be part of the simulation state:

* ``state.ud`` – displacement field with shape ``mesh.n + (3,)``.
* ``state.pd`` – momentum density field with shape ``mesh.n + (3,)``.
* ``state.material["C"]`` – stiffness tensor (Voigt notation, ``(6, 6)`` per cell).
* ``state.material["rho"]`` – mass density.
* ``state.material["lambda_100"]`` and ``state.material["lambda_111"]`` –
  magnetostriction constants used by the magnetoelastic field.

When the displacement or strain should be imposed externally you can supply a
callable to ``MagnetoElasticField(ud=...)`` or ``MagnetoElasticField(
mechanical_strain=...)`` instead of using ``state.ud``.

Boundary Conditions
===================

Boundary conditions for the elastic system are attached to the state through
``state.bcs``:

.. code-block:: python

    dirichlet = [BC(mask, lambda state: prescribed_ud)]
    traction = [PlaneBC(Plane(0, 0, -1), lambda state: surface_force)]
    state.bcs = {"ud": dirichlet, "t": traction}

``BC`` fixes the displacement on selected grid points, while
:class:`.PlaneBC` applies stresses on planar patches.  Several conditions can
be combined in the corresponding lists.

Self-Consistent Time Integration
================================

The :class:`.LLGWithLESolver` integrates magnetization dynamics together with
the elastic wave equation.  Provide the magnetoelastic term alongside the
usual micromagnetic fields:

.. code-block:: python

    mesh = Mesh(n, dx, pbc=(0, 1, 1))
    state = State(mesh)
    state.material = {
        "Ms": state.Constant(490e3),
        "A": state.Constant(8e-12),
        "rho": state.Constant(8900.0),
        "C": state.Constant(C_cubic(C11, C12, C44)),
        "lambda_100": state.Constant(-4.6e-5),
        "lambda_111": state.Constant(-2.4e-5),
        "alpha": state.Constant(0.05),
    }
    state.m = state.Constant((0, 0, 1))
    state.ud = state.Constant((0, 0, 0))
    state.pd = state.Constant((0, 0, 0))

    magnetoelastic = MagnetoElasticField()
    exchange = ExchangeField()
    external = ExternalField(h0 * state.Constant((0, 0, 1)))

    llg = LLGWithLESolver(
        [magnetoelastic, exchange, external],
        C_sym="cubic",
        boundary_nodes=3,
    )

The ``C_sym`` argument narrows the stiffness tensor to a specific symmetry
class which reduces the number of elastic coefficients that have to be
evaluated.  The ``boundary_nodes`` parameter selects the stencil that is used
to impose traction boundary conditions.

Example: Driven Shear Wave
==========================

The demo ``demos/linear_elasticity/bulk_shear_mode`` drives a shear wave
into a permalloy stripe by prescribing a sinusoidal traction on one end and a
fixed displacement on the opposite boundary.  The magnetoelastic field turns
the elastic energy into an effective magnetic field which, in combination
with precession and damping, leads to coupled spin–lattice oscillations.

As in the example, it is convenient to log the displacement, momentum density
and magnetization simultaneously:

.. code-block:: python

    logger = Logger(
        "log",
        scalars=["t", llg.E, llg.U, llg.T_el],
        fields=["m", "ud", "pd"],
        fields_every=N_T,
    )
    for _ in range(num_steps):
        logger << state
        llg.step(state, dt)

Complete Code
=============

The full script used for the shear-wave experiment can be found at
:download:`bulk_shear_mode/run.py <../demos/linear_elasticity/bulk_shear_mode/run.py>`.

Reference
=========

.. autoclass:: LLGWithLESolver
