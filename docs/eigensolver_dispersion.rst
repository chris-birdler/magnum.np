.. currentmodule:: magnumnp

##############################################
Spin Wave Dispersion (run_dispersion.py)
##############################################

This example uses the eigensolver to calculate the spin wave dispersion relation in a magnonic waveguide, following [VENKAT]_.

Setup
=====

The waveguide is discretized on a 400 x 20 x 1 grid with a cell size of 2.5 nm x 2.5 nm x 1 nm.
The material parameters are :math:`M_s = 800\,\text{kA/m}` and :math:`A = 13\,\text{pJ/m}`.
An external bias field of 804 kA/m is applied along the x-direction:

.. code-block:: python

    n = (400, 20, 1)
    dx = (2.5e-9, 2.5e-9, 1e-9)

    mesh = Mesh(n, dx)
    state = State(mesh)
    state.material = {
            "A": 13e-12,
            "Ms": 800e3,
            }

    state.m = state.Constant([1.,0.,0.])

    demag    = DemagField()
    exchange = ExchangeField()
    external = ExternalField([804e3, 0.0, 0.0])

Groundstate
===========

The equilibrium magnetization is obtained by energy minimization:

.. code-block:: python

    minimizer = MinimizerBB([demag, exchange, external])
    minimizer.minimize(state)

Eigenmode Calculation
=====================

A large number of eigenmodes (200) is computed to resolve the dispersion branches:

.. code-block:: python

    eigen = EigenSolver(state, [demag, exchange], [external])
    res = eigen.solve(k=200, tol=1e-6)

Dispersion Relation
===================

Since the eigensolver already works in the frequency domain, only the spatial eigenvectors need to be Fourier transformed to obtain the dispersion relation :math:`\omega(k)`.
The user specifies the line along which the spatial FFT is performed and the magnetization component to use, via a ``lambda`` function:

.. code-block:: python

    points = lambda vvv: vvv[:,10,0,2,:]
    omega = res.omega.cpu().numpy()
    kk, ww, mm = res.dispersion(points, dx[0])

This selects the :math:`m_z` component along the center of the waveguide (:math:`y = 10`).
The result is plotted as a 2D color map of the spectral weight in :math:`(k_x, f)` space:

.. code-block:: python

    fig, ax = plt.subplots()
    ax.pcolormesh(kk*1e-9, ww*1e-9/2./np.pi, np.log10(np.abs(mm)**2.), cmap="viridis")
    ax.set_xlim([-0.3, 0.3])
    ax.set_ylim([40, 80])
    ax.set_xlabel('$k_{x}\,\mathrm{[rad/nm]}$')
    ax.set_ylabel('Frequency f [GHz]')

.. image:: _static/dispersion.png
   :width: 400


Complete Code
=============

The complete code can be viewed here: :download:`run_dispersion.py <../demos/eigensolver/run_dispersion.py>`.


References
==========

.. [VENKAT] G. Venkat, H. Fangohr, and A. Prabhakar, "Absorbing boundary layers for spin wave micromagnetics," Journal of Magnetism and Magnetic Materials 450, 34--39 (2018).
