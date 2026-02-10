.. currentmodule:: magnumnp

###########
Eigensolver
###########

The eigenmode solver of *magnum.np* computes spin-wave resonances and their spectral signatures directly from a linearization of the Landau--Lifshitz--Gilbert (LLG) equation.
The undamped eigenmodes are obtained from a Hermitian generalized eigenvalue problem solved with the Lanczos algorithm.
Damping is incorporated via first-order perturbation theory, and closed-form expressions for the magnetization power spectrum and the absorbed magnetic power are derived using a modal expansion.
All continuum formulas are written in terms of a volume-averaged inner product, so that they translate directly to the finite-difference discretization without rescaling.

.. toctree::
   :maxdepth: 1

   Elliptical Nanodisc <eigensolver_carlotti>
   AFM Coupled Bilayer <eigensolver_bilayer>
   Spin Wave Dispersion <eigensolver_dispersion>

Eigenmode Equation
==================

The following calculations are based on the work of d'Aquino et al. [dAquino2009]_.
Linearization of the LLG equation for small fluctuations :math:`\boldsymbol{v}(\boldsymbol{x},t)` around a stable equilibrium :math:`\boldsymbol{m}_0` leads to

.. math::

   -\alpha\,\dot{\boldsymbol{v}} - \boldsymbol{m}_0 \times \dot{\boldsymbol{v}}
   = \underbrace{(\mathbb{1} - \boldsymbol{m}_0 \otimes \boldsymbol{m}_0)}_{\mathcal{P}_0}\;
     \underbrace{\gamma\!\left(h_0\,\mathbb{1} - \boldsymbol{h}^{\mathrm{lin}}\right)}_{\mathcal{A}_0}[\boldsymbol{v}],

with the gyromagnetic ratio :math:`\gamma = 2.2127615 \cdot 10^5\,\text{m/As}`, the dimensionless damping parameter :math:`\alpha`, the tangent-plane projection operator :math:`\mathcal{P}_0`, the linearized effective field operator :math:`\boldsymbol{h}^{\mathrm{lin}}[\boldsymbol{v}] = \frac{\delta \boldsymbol{h}^{\mathrm{eff}}}{\delta \boldsymbol{m}}[\boldsymbol{v}]`, and the parallel component of the equilibrium field :math:`h_0 = \boldsymbol{m}_0 \cdot \boldsymbol{h}^{\mathrm{eff}}[\boldsymbol{m}_0]`.

Using the Ansatz :math:`\boldsymbol{v}(\boldsymbol{x},t) = \sum_k a_k\,\tilde{\boldsymbol{\varphi}}'_k(\boldsymbol{x})\,e^{i\omega'_k t}` with the damped eigenvectors :math:`\tilde{\boldsymbol{\varphi}}'_k` and eigenfrequencies :math:`\omega'_k`, the differential equation transforms into the algebraic generalized eigenvalue equation:

.. math::

   \omega'_k\underbrace{(-i\alpha\,\tilde{\boldsymbol{\varphi}}'_k)}_{i\,\delta\mathcal{B}[\tilde{\boldsymbol{\varphi}}'_k]}
   + \omega'_k\underbrace{(-i\boldsymbol{m}_0 \times \tilde{\boldsymbol{\varphi}}'_k)}_{\mathcal{B}_0[\tilde{\boldsymbol{\varphi}}'_k]}
   = \underbrace{\mathcal{P}_0\,\mathcal{A}_0[\tilde{\boldsymbol{\varphi}}'_k]}_{\mathcal{A}_{0\perp}[\tilde{\boldsymbol{\varphi}}'_k]}.


Without Damping (:math:`\alpha = 0`)
-------------------------------------

Setting :math:`\alpha = 0` yields the Hermitian generalized eigenvalue problem

.. math::

   \mathcal{A}_{0\perp}[\tilde{\boldsymbol{\varphi}}_k] = \omega_k\,\mathcal{B}_0[\tilde{\boldsymbol{\varphi}}_k],

where :math:`\tilde{\boldsymbol{\varphi}}_k` and :math:`\omega_k` denote the undamped eigenmodes and eigenfrequencies, respectively.
Both operators :math:`\mathcal{A}_{0\perp}` and :math:`\mathcal{B}_0` are Hermitian, which allows using the Lanczos algorithm for the solution of the eigenvalue problem, resulting in more stable and efficient calculations.

The volume-averaged inner product is defined as

.. math::

   \langle f,g \rangle = \frac{1}{V}\int f^* \cdot g\,\mathrm{d}\boldsymbol{x}, \qquad
   \langle f,g \rangle_{\mathcal{O}} = \langle f,\,\mathcal{O}\,g \rangle, \qquad
   \|f\|^2 = \langle f,f \rangle.

In the discrete counterpart implemented in *magnum.np* this becomes

.. math::

   \langle f,g \rangle = \frac{1}{N_{\mathrm{mag}}}\sum_{i\in\Omega_{\mathrm{mag}}} f_i^*\,g_i,

which is a mean over the magnetic domain.
Because all continuum formulas are expressed using :math:`\langle\cdot,\cdot\rangle`, the discrete expressions take exactly the same form without requiring additional prefactors.

The orthogonality relation reads

.. math::

   \langle \tilde{\boldsymbol{\varphi}}_h, \tilde{\boldsymbol{\varphi}}_k \rangle_{\mathcal{A}_{0\perp}} = \omega_h\;\langle \tilde{\boldsymbol{\varphi}}_h, \tilde{\boldsymbol{\varphi}}_k \rangle_{\mathcal{B}_0} = \delta_{hk}.


Perturbation Analysis for Damping (:math:`\alpha > 0`)
-------------------------------------------------------

For small damping :math:`\alpha > 0`, a first-order perturbation analysis with :math:`\delta\mathcal{B} = -\alpha\,\mathbb{1}` yields the frequency correction

.. math::

   \delta\omega_k = -\omega_k\,\frac{\langle \tilde{\boldsymbol{\varphi}}_k,\tilde{\boldsymbol{\varphi}}_k \rangle_{\delta\mathcal{B}}}{\langle \tilde{\boldsymbol{\varphi}}_k,\tilde{\boldsymbol{\varphi}}_k \rangle_{\mathcal{B}_0}} = \omega_k^2\,\langle \tilde{\boldsymbol{\varphi}}_k,\alpha\,\tilde{\boldsymbol{\varphi}}_k \rangle.

The perturbation analysis results in an additional imaginary contribution to the eigenfrequency, which leads to a damped harmonic oscillation in the time domain.
Ignoring the small perturbations of the eigenvectors :math:`\tilde{\boldsymbol{\varphi}}_k` one ends up with

.. math::

   \boldsymbol{v}_\alpha(\boldsymbol{x},t) \approx \sum_k a_k\,\tilde{\boldsymbol{\varphi}}_k(\boldsymbol{x})\,e^{i\omega_k t}\,e^{-\delta\omega_k t}.


Harmonic Excitation
-------------------

Adding :math:`\tilde{\boldsymbol{h}}^{\mathrm{a}}` as a source term to the unperturbed eigenvalue problem yields

.. math::

   \omega\,\mathcal{B}_0\,\delta\hat{\boldsymbol{m}} = \mathcal{A}_{0\perp}\,\delta\hat{\boldsymbol{m}} - \gamma\,\mathcal{P}_0\,\tilde{\boldsymbol{h}}^{\mathrm{a}}.

The excited solution :math:`\delta\hat{\boldsymbol{m}}` is expressed by the unperturbed eigenvectors

.. math::

   \delta\hat{\boldsymbol{m}}(\boldsymbol{x},\omega) = \sum_k \tilde{a}_k(\omega)\,\tilde{\boldsymbol{\varphi}}_k(\boldsymbol{x}).

Using the orthogonality relation to calculate the mode amplitudes gives

.. math::

   \tilde{a}_h = \frac{\gamma\,\omega_h}{\omega'_h - \omega}\,\tilde{h}_h, \qquad
   \tilde{h}_h = \langle \tilde{\boldsymbol{\varphi}}_h, \mathcal{P}_0\,\tilde{\boldsymbol{h}}^{\mathrm{a}} \rangle,

where the damped eigenfrequencies :math:`\omega'_h = \omega_h + i\,\delta\omega_h` are used in the denominator.


Power Spectrum of Magnetization
--------------------------------

The volume-averaged power spectrum of the magnetization driven by a harmonic excitation :math:`\tilde{\boldsymbol{h}}^{\mathrm{a}}` is defined as (compare Eqn. (22) of Ref. [dAquino2023]_)

.. math::

   \tilde{p}(\omega) = \frac{1}{2V}\int M_s^2\,|\delta\hat{\boldsymbol{m}}(\boldsymbol{x},\omega)|^2\;\mathrm{d}\boldsymbol{x}
   = \frac{1}{2}\|M_s\,\delta\hat{\boldsymbol{m}}(\boldsymbol{x},\omega)\|^2,

with :math:`\delta\hat{\boldsymbol{m}}` being the harmonic magnetization response driven by :math:`\tilde{\boldsymbol{h}}^{\mathrm{a}}`.
Inserting the modal expansion yields

.. math::

   \tilde{p}(\omega) = \frac{\gamma^2}{2}\sum_{k,h}\frac{\omega_k\,\omega_h}{(\omega'_k - \omega)(\omega'^*_h - \omega)}\;\tilde{h}_k\,\tilde{h}_h^*\,\langle M_s\,\tilde{\boldsymbol{\varphi}}_k,\,M_s\,\tilde{\boldsymbol{\varphi}}_h \rangle.

For well separated resonances the off-diagonal terms (:math:`k \neq h`) can be neglected, giving

.. math::

   \tilde{p}(\omega) \approx \frac{\gamma^2}{2}\sum_k \frac{\omega_k^2\,|\tilde{h}_k|^2\,\|M_s\,\tilde{\boldsymbol{\varphi}}_k\|^2}{(\omega_k - \omega)^2 + \delta\omega_k^2}
   = \frac{1}{2}\sum_k |\tilde{a}_k|^2\,\|M_s\,\tilde{\boldsymbol{\varphi}}_k\|^2.


Absorbed Magnetic Power
------------------------

The average power absorbed by the magnetic system :math:`P_{\mathrm{abs}}(\omega)` driven by a harmonic excitation :math:`\tilde{\boldsymbol{h}}^{\mathrm{a}}` can be defined as the real part of the complex absorbed power :math:`\hat{P}_{\mathrm{abs}}(\omega)` (compare Eq. (25) of Ref. [dAquino2023]_)

.. math::

   \hat{P}_{\mathrm{abs}}(\omega) = \frac{\mu_0}{2V}\int i\omega\,M_s(\boldsymbol{x})\,\delta\hat{\boldsymbol{m}}(\boldsymbol{x},\omega) \cdot \tilde{\boldsymbol{h}}^{\mathrm{a}*}(\boldsymbol{x})\;\mathrm{d}\boldsymbol{x}
   = \frac{\mu_0}{2}\,i\omega\,\langle \tilde{\boldsymbol{h}}^{\mathrm{a}},\,M_s\,\delta\hat{\boldsymbol{m}} \rangle.

Inserting the modal expansion and defining the :math:`M_s`-weighted modal projection :math:`\tilde{g}_k = \langle \tilde{\boldsymbol{\varphi}}_k,\,M_s\,\mathcal{P}_0\,\tilde{\boldsymbol{h}}^{\mathrm{a}} \rangle` yields

.. math::

   \hat{P}_{\mathrm{abs}}(\omega) = \frac{\mu_0}{2}\sum_k i\omega\,\tilde{a}_k\,\tilde{g}_k^*
   = \frac{\mu_0\,\gamma}{2}\sum_k \frac{i\omega\,\omega_k\;\tilde{h}_k\,\tilde{g}_k^*}{\omega_k - \omega + \delta\omega_k}.


References
==========

.. [dAquino2009] M. d'Aquino, C. Serpico, G. Miano, and C. Forestiere, "A novel formulation for the numerical computation of magnetization modes in complex micromagnetic systems," Journal of Computational Physics 228, 6130--6149 (2009).

.. [dAquino2023] M. d'Aquino and R. Hertel, "Micromagnetic frequency-domain simulation methods for magnonic systems," Journal of Applied Physics 133, 033902 (2023).
