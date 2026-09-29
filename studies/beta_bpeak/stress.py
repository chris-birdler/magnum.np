"""
Random internal stress field (model of dislocation stress fields).

Model (assumption, not derived from micromechanics):
  * each of the 6 stress components (Voigt order xx, yy, zz, yz, xz, xy)
    is an independent, stationary Gaussian random field,
  * zero mean, standard deviation sigma_rms (over the magnetic cells),
  * Gaussian covariance  C(r) = sigma_rms^2 * exp(-r^2 / (2 xi^2)),
  * periodic (it is made with FFTs on the periodic box).

The field does NOT satisfy mechanical equilibrium (div sigma = 0). The
hydrostatic part has no magnetic effect (isotropic magnetostriction), see
field_terms_extra.StressAnisotropyField.

Implementation: white noise w(r) is filtered with the kernel whose Fourier
transform is G(k) = exp(-k^2 xi^2 / 4). The filtered field then has the
covariance exp(-r^2 / (2 xi^2)) (up to a constant, removed by the
normalisation to sigma_rms).
"""
import math
import warnings
import numpy as np

__all__ = ["random_stress_field"]


def random_stress_field(N, dx, xi, sigma_rms, seed, mask=None):
    """
    :param N:         cells per edge
    :param dx:        cell size [m]
    :param xi:        correlation length [m]
    :param sigma_rms: standard deviation of each component [Pa]
    :param seed:      RNG seed
    :param mask:      bool array [N,N,N]; the statistics are normalised on
                      these cells (the magnetic cells). None = all cells.
    :returns: float64 array [N,N,N,6] in Pa (Voigt order xx,yy,zz,yz,xz,xy)
    """
    if sigma_rms == 0.0:
        return np.zeros((N, N, N, 6))
    if xi < 3.0 * dx:
        warnings.warn("xi = %.3g m is less than 3 cells (dx = %.3g m): the stress field "
                      "is not resolved and grid pinning can dominate." % (xi, dx))

    rng = np.random.default_rng(seed)
    k = 2.0 * math.pi * np.fft.fftfreq(N, d=dx)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    G = np.exp(-(kx**2 + ky**2 + kz**2) * xi**2 / 4.0)

    if mask is None:
        mask = np.ones((N, N, N), dtype=bool)

    sig = np.empty((N, N, N, 6))
    for c in range(6):
        w = rng.standard_normal((N, N, N))
        f = np.fft.ifftn(np.fft.fftn(w) * G).real
        f -= f[mask].mean()
        f *= sigma_rms / f[mask].std()
        sig[..., c] = f
    return sig
