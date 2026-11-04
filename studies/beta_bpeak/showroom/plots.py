"""Showroom loss equation (Chris 2026-11-04): trend plots of
P = [1 + p_c (kappa_0 + kappa_a a_c)] [c_m f K_eff (1 + gamma_h h) (1 + gamma_s phi_s) (B/J_s)^beta [1 + (f/f_r)^(alpha-1)]
                                      + pi^2 d^2 f^2 B^2 / (20 rho)]
(f_r = gamma J_s / (3 pi (mu_r - 1)) the Snoek resonance with mu_r - 1 = J_s^2 / (3 mu_0 K_eff) of the same
random-anisotropy rotation model, thus f_r = gamma mu_0 K_eff / (pi J_s); phi_s = 1 - (1 - 2 delta_s/d)^3 the volume
fraction of the surface layer; h = dK/K_eff the inhomogeneity)
One input varies, the others stay at the base point. y = P / P_base. Linear axes, no base-point mark, only a real
minimum of the total loss is marked (Chris 2026-11-04). The same base point and ranges as the deck in slide_decks,
folder showroom_loss_equation.
Writes fig_<name>_<theme>.png for theme light and dark."""
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent

GAMMA = 1.76e11          # gyromagnetic ratio (rad/(s T))
MU0 = 4e-7 * math.pi
# Base point for the MHz comparison (Chris 2026-11-04): f 30 MHz in the range 1 to 50 MHz, B 10 mT, particle material
# with mu_r about 30 (K_eff 20 kJ/m^3 with J_s 1.5 T: f_r 939 MHz), iron-nickel-like (rho 0.43 uOhm m, Kanthal
# Nifethal 52), d 3 um (eddy current part 30 % of the loss). c_m 0.048 and gamma_s 0.26 give 4.4 W/cm^3 at the base
# point (ten times 0.44 W/cm^3, Chris: from measurements) and d_opt 3.3 um at 1 MHz (Chris: 3 to 4 um).
BASE = dict(f=30e6, alpha=1.5, B=0.01, Js=1.5, K=20e3, h=0.2, gh=2.0, ds=0.1e-6, gs=0.26, d=3e-6, rho=0.43e-6,
            beta=2.0, c_h=0.048, p_c=0.2, a_c=0.05, k0=1.0, ka=20.0)

THEMES = {
    "light": dict(fg="#1d2733", muted="#5b6573", grid="#d9dee5", hyst="#2a78d6", eddy="#eb6834", total="#1d2733",
                  ramp=["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]),
    "dark": dict(fg="#eef1f4", muted="#aab3bf", grid="#36404c", hyst="#3987e5", eddy="#d95926", total="#eef1f4",
                 ramp=["#1c5cab", "#3987e5", "#86b6ef", "#cde2fb"]),
}


def phi_s(p):
    """Volume fraction of the surface layer of thickness delta_s in a sphere of diameter d (1 if 2 delta_s >= d)."""
    x = 2.0 * p["ds"] / p["d"]
    return 1.0 - (1.0 - x) ** 3 if x < 1.0 else 1.0


def mu_r(p):
    """Relative permeability of the particle material from K_eff and J_s (rotation with random easy axes)."""
    return 1.0 + p["Js"] ** 2 / (3.0 * MU0 * p["K"])


def f_r(p):
    """Snoek resonance frequency (Hz) with mu_r from K_eff: gamma mu_0 K_eff / (pi J_s)."""
    return GAMMA * p["Js"] / (3.0 * math.pi * (mu_r(p) - 1.0))


def skin_depth(p):
    """Skin depth (m) in the particle material with mu_r below f_r."""
    return math.sqrt(2.0 * p["rho"] / (2.0 * math.pi * p["f"] * MU0 * mu_r(p)))


def parts(p):
    """(magnetic part, eddy current part, contact factor) in W/m^3 and 1."""
    ph = (p["c_h"] * p["f"] * p["K"] * (1.0 + p["gh"] * p["h"]) * (1.0 + p["gs"] * phi_s(p))
          * (p["B"] / p["Js"]) ** p["beta"] * (1.0 + (p["f"] / f_r(p)) ** (p["alpha"] - 1.0)))
    pe = math.pi ** 2 * p["d"] ** 2 * p["f"] ** 2 * p["B"] ** 2 / (20.0 * p["rho"])
    fc = 1.0 + p["p_c"] * (p["k0"] + p["ka"] * p["a_c"])
    return ph, pe, fc


P0 = (lambda h, e, c: c * (h + e))(*parts(BASE))


def d_opt(p):
    """Particle diameter of the smallest loss (numerical minimum of the magnetic and the eddy current part over d)."""
    ds = np.geomspace(1e-8, 1e-3, 20001)
    tot = [sum(parts(dict(p, d=x))[:2]) for x in ds]
    return float(ds[int(np.argmin(tot))])


def sweep(key, xs, **over):
    h, e, t = [], [], []
    for x in xs:
        p = dict(BASE, **over)
        p[key] = x
        ph, pe, fc = parts(p)
        h.append(fc * ph / P0)
        e.append(fc * pe / P0)
        t.append(fc * (ph + pe) / P0)
    return np.array(h), np.array(e), np.array(t)


# name, key, range (SI), axis scale factor, x label, x ticks (axis units), label of a marked minimum.
# K_eff 5 to 60 kJ/m^3 keeps f_r above the base frequency; f up to 50 MHz stays below f_r; d up to 5 um stays below
# the skin depth at 30 MHz (11 um).
PANELS = [
    ("K", "K", np.linspace(5e3, 60e3, 221), 1e-3, r"effective anisotropy $K_\mathrm{eff}$ (kJ/m$^3$)",
     (5, 10, 20, 30, 40, 50, 60), None),
    ("Js", "Js", np.linspace(1.0, 2.0, 201), 1.0, r"saturation polarization $J_\mathrm{s}$ (T)",
     (1.0, 1.2, 1.4, 1.6, 1.8, 2.0), None),
    ("B", "B", np.linspace(0.002, 0.04, 191), 1e3, r"peak flux density $B_\mathrm{peak}$ (mT)",
     (2, 10, 20, 30, 40), None),
    ("f", "f", np.linspace(1e6, 50e6, 246), 1e-6, r"frequency $f$ (MHz)", (1, 10, 20, 30, 40, 50), None),
    ("d", "d", np.linspace(0.3e-6, 5e-6, 236), 1e6, r"particle diameter $d$ ($\mu$m)", (0.3, 1, 2, 3, 4, 5),
     r" $d_\mathrm{opt}$ %.1f $\mu$m"),
    ("rho", "rho", np.linspace(0.1e-6, 2e-6, 191), 1e6, r"resistivity $\rho$ ($\mu\Omega\,$m)",
     (0.1, 0.5, 1.0, 1.5, 2.0), None),
]


def style(ax, t):
    ax.set_facecolor("none")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["muted"])
    ax.tick_params(which="both", colors=t["muted"], labelsize=9)
    ax.grid(True, which="major", color=t["grid"], lw=0.7)
    ax.xaxis.label.set_color(t["fg"])
    ax.yaxis.label.set_color(t["fg"])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))


def new_fig():
    fig, ax = plt.subplots(figsize=(4.6, 3.2), dpi=200)
    fig.patch.set_alpha(0.0)
    return fig, ax


def save(fig, name, theme):
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / ("fig_%s_%s.png" % (name, theme)), transparent=True)
    plt.close(fig)


def legend(ax, t, loc):
    lg = ax.legend(frameon=False, fontsize=8.5, loc=loc)
    for tx in lg.get_texts():
        tx.set_color(t["fg"])


def parts_lines(ax, t, x, h, e, tot):
    """Magnetic and eddy current part, and the total dashed on top so that no part hides under it."""
    ax.plot(x, h, color=t["hyst"], lw=1.6, label="magnetic part")
    ax.plot(x, e, color=t["eddy"], lw=1.6, label="eddy current part")
    ax.plot(x, tot, color=t["total"], lw=2.0, ls=(0, (4, 2)), zorder=3, label="total")


def power_panels(theme):
    t = THEMES[theme]
    for name, key, xs, sc, xl, xt, min_label in PANELS:
        h, e, tot = sweep(key, xs)
        x = xs * sc
        fig, ax = new_fig()
        style(ax, t)
        parts_lines(ax, t, x, h, e, tot)
        i = int(np.argmin(tot))
        if min_label and 0 < i < len(tot) - 1:                 # only a real minimum, not an end of the range
            ax.axvline(x[i], color=t["muted"], lw=0.8, ls=":")
            ax.plot([x[i]], [tot[i]], "o", ms=5, color=t["total"], zorder=4)
            ax.text(x[i], 0.04, min_label % x[i], transform=ax.get_xaxis_transform(), color=t["muted"],
                    fontsize=9, ha="left", va="bottom")
        ax.set_xlim(x[0], x[-1])
        ax.set_xticks(xt)
        ax.set_ylim(0, 1.05 * tot.max())
        ax.set_xlabel(xl, fontsize=10)
        ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
        legend(ax, t, "best")
        save(fig, name, theme)


def homogeneity_panel(theme):
    t = THEMES[theme]
    fig, ax = new_fig()
    style(ax, t)
    hs = np.linspace(0.0, 1.0, 101)
    h, e, tot = sweep("h", hs)
    parts_lines(ax, t, hs, h, e, tot)
    ax.set_xlabel(r"inhomogeneity $h$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05 * tot.max())
    legend(ax, t, "upper left")
    save(fig, "h", theme)


def contact_panels(theme):
    t = THEMES[theme]
    # loss against the contact probability, one line per contact area fraction
    fig, ax = new_fig()
    style(ax, t)
    ps = np.linspace(0.0, 1.0, 101)
    for col, a in zip(t["ramp"], (0.0, 0.05, 0.1, 0.2)):
        _, _, tot = sweep("p_c", ps, a_c=a)
        ax.plot(ps, tot, color=col, lw=1.8, label=r"$a_\mathrm{c}$ %.2f" % a)
    ax.set_xlabel(r"contact probability $p_\mathrm{c}$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    legend(ax, t, "upper left")
    save(fig, "pc", theme)
    # loss against the contact area fraction, one line per contact probability
    fig, ax = new_fig()
    style(ax, t)
    acs = np.linspace(0.0, 0.3, 101)
    for col, pc in zip(t["ramp"], (0.0, 0.2, 0.5, 1.0)):
        _, _, tot = sweep("a_c", acs, p_c=pc)
        ax.plot(acs, tot, color=col, lw=1.8, label=r"$p_\mathrm{c}$ %.1f" % pc)
    ax.set_xlabel(r"contact area fraction $a_\mathrm{c}$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 0.3)
    ax.set_ylim(0, None)
    legend(ax, t, "upper left")
    save(fig, "ac", theme)


if __name__ == "__main__":
    plt.rcParams["mathtext.fontset"] = "cm"
    for th in THEMES:
        power_panels(th)
        homogeneity_panel(th)
        contact_panels(th)
    ph, pe, fc = parts(BASE)
    print("base: magnetic %.3g W/m3, eddy %.3g W/m3, contact factor %.3f, P_base %.3g W/m3 = %.3g W/cm3"
          % (ph, pe, fc, P0, P0 * 1e-6))
    print("mu_r %.3g; f_r %.0f MHz; share magnetic %.2f, eddy %.2f; phi_s %.4f; d_opt %.2f um at 30 MHz, %.2f um at "
          "1 MHz; skin depth %.1f um" % (mu_r(BASE), f_r(BASE) / 1e6, ph / (ph + pe), pe / (ph + pe), phi_s(BASE),
                                         1e6 * d_opt(BASE), 1e6 * d_opt(dict(BASE, f=1e6)), 1e6 * skin_depth(BASE)))
