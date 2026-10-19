"""Standalone evaluation of per-powder BH-loop logs (numpy only).

Evaluates a .dat file produced by the logger set up in per_powder_logging.py
and computes, for each powder phase over one closed field cycle:

    B_peak_x  peak flux density in the drive direction, using the FULL
              local field (self-demag included -- flux continuity requires
              it):  B_i = mu0 * (H_ext + <h_demag>_i) + J_i
    J_peak    peak |J_i|; sanity check, must be <= mu0 * Ms_i
    w_ext     loop integral of H_ext . dJ_i
              -> work put into phase i by the external source. The
                 phase-volume-weighted sum over both phases equals the
                 total composite loss EXACTLY, but it ignores
                 magnetostatic energy transfer between the powders.
    w_loc     loop integral of (H_ext + <h_demag>_i) . dJ_i
              -> approximates the energy actually DISSIPATED in phase i.
                 The self-demag contribution cancels over a closed cycle
                 (the self-energy is a state function); the residual error
                 stems from using products of phase averages instead of
                 averages of products.
    tr_in     loop integral of <h_cross>_i . dJ_i
              -> net magnetostatic energy pumped into phase i by the
                 OTHER powder per cycle.
    w_self    loop integral of <h_self>_i . dJ_i
              -> should be ~ 0 (conservative term); its magnitude shows
                 the size of the phase-averaging error.

Consistency checks printed at the end:
    * phi_L*w_ext_L + phi_S*w_ext_S must equal the total composite loss
      (exact identity, independent of any approximation).
    * phi_L*tr_L + phi_S*tr_S ~ 0 (energy leaving one powder enters the
      other).

All loss numbers are per unit PHASE volume and per cycle [J/m^3]; multiply
by phi_i * V_sample for absolute energies. The evaluation window must span
exactly one full field period (steady state, i.e. exclude the initial
transient), otherwise the conservative terms do not cancel.

Usage::

    python evaluate_bh_loop.py data/loop.dat --period 5e-9
    python evaluate_bh_loop.py data/loop.dat --period 5e-9 --t-end 10e-9
"""

import numpy as np

MU0 = 4e-7 * np.pi  # vacuum permeability [Vs/Am]; identical to magnumnp.constants.mu_0

# np.trapz was renamed to np.trapezoid in NumPy 2.0; support both
_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))


def read_log(path):
    """Parse a ScalarLogger .dat file into a {column_name: 1-D array} dict."""
    with open(path) as f:
        header = f.readline().lstrip('#').split()
    data = np.loadtxt(path)
    if data.ndim == 1:                      # single-row file
        data = data.reshape(1, -1)
    return {name: data[:, i] for i, name in enumerate(header)}


def _vec(cols, name):
    """Stack the _x/_y/_z columns of a vector quantity into shape [N, 3]."""
    return np.stack([cols[name + '_' + c] for c in 'xyz'], axis=1)


def loop_integral(h, j):
    """Closed-cycle integral  w = loop integral of h . dJ  [J/m^3 per cycle].

    h [A/m] and j [T] are time series of shape [N, 3] sampled along the
    loop; the trapezoidal rule integrates each vector component along the
    path and sums them (with the drive along x the transverse components
    are usually negligible, but including them costs nothing).

    The window MUST span exactly one full field period (first and last
    sample at the same field phase), otherwise the conservative demag
    self-term does not cancel and contaminates the result.
    """
    return sum(_trapz(h[:, c], j[:, c]) for c in range(3))


def evaluate(path, period, t_end=None, verbose=True):
    """Evaluate per-powder B_peak, losses and transfer over one field period.

    :param path:    ScalarLogger .dat file written by per_powder_logging.py
    :param period:  duration of one full field cycle [s]
    :param t_end:   end of the evaluation window [s]; defaults to the last
                    logged time. Choose the LAST period so that the initial
                    transient (first magnetisation curve) is excluded.
    :return:        nested dict with per-phase results and global checks
    """
    cols = read_log(path)
    t = cols['t']
    if t_end is None:
        t_end = t[-1]
    sel = (t >= t_end - period - 1e-15) & (t <= t_end + 1e-15)
    if sel.sum() < 16:
        raise ValueError("Fewer than 16 samples in the selected period -- "
                         "log more densely or check the --period value.")

    h_ext = _vec(cols, 'h_ext')[sel]
    phi = {'L': cols['phi_L'][0], 'S': cols['phi_S'][0]}
    cross = {'L': 'hS_L', 'S': 'hL_S'}      # field of the OTHER powder in phase i
    own = {'L': 'hL_L', 'S': 'hS_S'}        # field of phase i in itself

    res, J_tot = {}, 0.0
    for ph in ('L', 'S'):
        J  = _vec(cols, 'J_' + ph)[sel]
        hd = _vec(cols, 'hd_' + ph)[sel]
        h_cross = _vec(cols, cross[ph])[sel]
        h_own   = _vec(cols, own[ph])[sel]

        # flux density with the FULL local field (self-demag included)
        B = MU0 * (h_ext + hd) + J

        res[ph] = {
            'B_peak_x': np.max(np.abs(B[:, 0])),            # [T], drive direction
            'J_peak':   np.max(np.linalg.norm(J, axis=1)),  # [T], must be <= mu0*Ms_i
            'w_ext':    loop_integral(h_ext, J),            # source work    [J/m^3 phase]
            'w_loc':    loop_integral(h_ext + hd, J),       # ~ dissipation  [J/m^3 phase]
            'tr_in':    loop_integral(h_cross, J),          # transfer from other powder
            'w_self':   loop_integral(h_own, J),            # ~ 0 check (conservative term)
        }
        J_tot = J_tot + phi[ph] * J        # composite polarisation (matrix has M=0)

    # --- global consistency checks -----------------------------------------
    w_tot = loop_integral(h_ext, J_tot)    # exact total loss [J/m^3 TOTAL volume]
    checks = {
        # must match w_tot exactly (work decomposition is complete):
        'sum_w_ext':        phi['L'] * res['L']['w_ext'] + phi['S'] * res['S']['w_ext'],
        'w_tot':            w_tot,
        # must be ~ 0 (energy leaving one powder enters the other):
        'transfer_balance': phi['L'] * res['L']['tr_in'] + phi['S'] * res['S']['tr_in'],
        # each must be ~ 0 (self-demag is conservative over a closed cycle;
        # residual = phase-averaging error):
        'w_self_L':         res['L']['w_self'],
        'w_self_S':         res['S']['w_self'],
    }

    if verbose:
        print(f"evaluation window: t = {t_end - period:.4e} .. {t_end:.4e} s "
              f"({int(sel.sum())} samples)")
        print(f"volume fractions:  phi_L = {phi['L']:.4f}, phi_S = {phi['S']:.4f}\n")
        hdr = f"{'quantity':<38}{'large (L)':>16}{'small (S)':>16}"
        print(hdr); print('-' * len(hdr))
        rows = [
            ('B_peak_x                       [T]', 'B_peak_x', '%16.4f'),
            ('J_peak  (sanity: <= mu0*Ms)    [T]', 'J_peak',   '%16.4f'),
            ('w_ext   (source work)     [J/m^3]',  'w_ext',    '%16.6e'),
            ('w_loc   (~ dissipation)   [J/m^3]',  'w_loc',    '%16.6e'),
            ('tr_in   (from other pow.) [J/m^3]',  'tr_in',    '%16.6e'),
            ('w_self  (should be ~ 0)   [J/m^3]',  'w_self',   '%16.6e'),
        ]
        for label, key, fmt in rows:
            print(f"{label:<38}" + fmt % res['L'][key] + fmt % res['S'][key])
        print("\nconsistency checks (per unit TOTAL volume):")
        print(f"  total loss  loop-int H_ext.dJ_tot  = {checks['w_tot']:.6e} J/m^3")
        print(f"  phi_L*w_ext_L + phi_S*w_ext_S      = {checks['sum_w_ext']:.6e} J/m^3  (must match)")
        print(f"  transfer balance (should be ~ 0)   = {checks['transfer_balance']:.6e} J/m^3")

    return {'phases': res, 'phi': phi, 'checks': checks}


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
        description="Evaluate per-powder B_peak and losses from a BH-loop log "
                    "written by the logger set up in per_powder_logging.py.")
    p.add_argument('logfile', help="path to the .dat file")
    p.add_argument('--period', type=float, required=True,
                   help="duration of one full field cycle in seconds")
    p.add_argument('--t-end', type=float, default=None,
                   help="end of the evaluation window in seconds "
                        "(default: last logged time)")
    a = p.parse_args()
    evaluate(a.logfile, a.period, t_end=a.t_end)
