"""
Reduced units of the study.

The physics of the model (T = 0, quasi-static LLG, periodic powder) depends
only on dimensionless numbers. Js and A only select the SI unit system:

    Ms   = Js / mu0                      [A/m]      field unit
    K_d  = mu0 Ms^2 / 2                  [J/m^3]    energy density unit
    l_ex = sqrt(2 A / (mu0 Ms^2))        [m]        length unit
    f_M  = gamma Ms / (2 pi)             [Hz]       frequency unit (28 GHz/T * Js)

Reduced parameters used by run_loops.py:
    d/l_ex, dx/l_ex, xi/l_ex        lengths
    phi                             packing fraction
    Q     = K / K_d                 anisotropy (uniaxial or cubic K1)
    kappa = K_sigma,rms / K_d       strength of the random stress anisotropy
    f/f_M                           drive frequency
    h = H / Ms,  b = B / Js,  w = W / K_d   (per cycle)

Default unit system: Js = 1.5 T, A = 10 pJ/m  ->  l_ex = 3.34 nm,
K_d = 895 kJ/m^3, f_M = 42.0 GHz. The default f/f_M gives 30 MHz there.
"""
import math

MU0 = 1.2566370614e-6
GAMMA = 2.21276157e5          # magnum.np constants.gamma [m/(A s)]

JS_REF = 1.5                  # T
A_REF = 10e-12                # J/m
F_REF = 30e6                  # Hz


def units(Js=JS_REF, A=A_REF):
    Ms = Js / MU0
    Kd = 0.5 * MU0 * Ms**2
    return {"Js": Js, "A": A, "Ms": Ms, "Kd": Kd,
            "l_ex": math.sqrt(2.0 * A / (MU0 * Ms**2)),
            "f_M": GAMMA * Ms / (2.0 * math.pi),
            "t_M": 1.0 / (GAMMA * Ms)}


F_REL_DEFAULT = F_REF / units()["f_M"]      # = 7.14e-4 -> 30 MHz at Js = 1.5 T
