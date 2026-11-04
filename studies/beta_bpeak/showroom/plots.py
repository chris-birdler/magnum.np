"""Showroom loss equation (Chris 2026-11-04): trend plots of
P = [1 + p_c (kappa_0 + kappa_a a_c)] [c_m f K_eff (1 + gamma_h h) (1 + gamma_s phi_s) (B/J_s)^beta [1 + (f/f_r)^(alpha-1)]
                                      + pi^2 d^2 f^2 B^2 / (20 rho)]
(f_r = gamma J_s / (3 pi (mu_r - 1)) the Snoek resonance with mu_r - 1 = J_s^2 / (3 mu_0 K_eff) of the same
random-anisotropy rotation model, thus f_r = gamma mu_0 K_eff / (pi J_s); phi_s = 1 - (1 - 2 delta_s/d)^3 the volume fraction of the surface layer; K_eff one input: for nanocrystalline materials after the Herzer averaging; h = dK/K_eff the inhomogeneity)
One input varies, the others stay at the base point. y = P / P_base (trends, not absolute losses).
Writes fig_<name>_<theme>.png for theme light and dark."""
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent

GAMMA = 1.76e11          # gyromagnetic ratio (rad/(s T))
MU0 = 4e-7 * math.pi
BASE = dict(f=1e5, alpha=1.5, B=0.1, Js=1.25, K=4.0, h=0.2, gh=2.0, ds=0.1e-6, gs=1.2, d=10e-6, rho=1.2e-6, beta=2.0, c_h=1.0,
            p_c=0.2, a_c=0.05, k0=1.0, ka=20.0)

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


def parts(p):
    """(magnetic part, eddy current part, contact factor) in W/m^3 and 1."""
    ph = (p["c_h"] * p["f"] * p["K"] * (1.0 + p["gh"] * p["h"]) * (1.0 + p["gs"] * phi_s(p))
          * (p["B"] / p["Js"]) ** p["beta"] * (1.0 + (p["f"] / f_r(p)) ** (p["alpha"] - 1.0)))
    pe = math.pi ** 2 * p["d"] ** 2 * p["f"] ** 2 * p["B"] ** 2 / (20.0 * p["rho"])
    fc = 1.0 + p["p_c"] * (p["k0"] + p["ka"] * p["a_c"])
    return ph, pe, fc


P0 = (lambda h, e, c: c * (h + e))(*parts(BASE))


def d_opt(p):
    """Particle diameter of the smallest loss (numerical minimum of hysteresis + eddy current over d)."""
    ds = np.geomspace(1e-8, 1e-3, 20001)
    tot = [sum(parts(dict(p, d=x))[:2]) for x in ds]
    return float(ds[int(np.argmin(tot))])


def d_opt_thin(p):
    """Thin-layer form (delta_s << d): d^3 = 60 gamma_s delta_s c_m K_eff (1 + gamma_h h) rho [1 + (f/f_r)^(alpha-1)]
    / (pi^2 J_s^2 f), exponent beta = 2."""
    return (60.0 * p["gs"] * p["ds"] * p["c_h"] * p["K"] * (1.0 + p["gh"] * p["h"]) * p["rho"]
            * (1.0 + (p["f"] / f_r(p)) ** (p["alpha"] - 1.0)) / (math.pi ** 2 * p["Js"] ** 2 * p["f"])) ** (1.0 / 3.0)


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


# name, key, range (SI), axis scale factor, x label, x ticks (axis units)
PANELS = [
    ("K", "K", np.geomspace(0.3, 300, 200), 1.0, r"effective anisotropy $K_\mathrm{eff}$ (J/m$^3$)",
     (0.3, 1, 3, 10, 30, 100, 300)),
    ("Js", "Js", np.geomspace(0.6, 2.0, 200), 1.0, r"saturation polarization $J_\mathrm{s}$ (T)",
     (0.6, 0.8, 1.0, 1.25, 1.6, 2.0)),
    ("B", "B", np.geomspace(0.01, 0.5, 200), 1e3, r"peak flux density $B_\mathrm{peak}$ (mT)",
     (10, 20, 50, 100, 200, 500)),
    ("f", "f", np.geomspace(1e3, 4.5e5, 200), 1e-3, r"frequency $f$ (kHz)", (1, 3, 10, 30, 100, 300)),
    ("d", "d", np.geomspace(0.1e-6, 100e-6, 300), 1e6, r"particle diameter $d$ ($\mu$m)",
     (0.1, 0.3, 1, 3, 10, 30, 100)),
    ("rho", "rho", np.geomspace(0.2e-6, 5e-6, 200), 1e6, r"resistivity $\rho$ ($\mu\Omega\,$m)",
     (0.2, 0.5, 1, 2, 5)),
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


def new_fig():
    fig, ax = plt.subplots(figsize=(4.6, 3.2), dpi=200)
    fig.patch.set_alpha(0.0)
    return fig, ax


def save(fig, name, theme):
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / ("fig_%s_%s.png" % (name, theme)), transparent=True)
    plt.close(fig)


def power_panels(theme):
    t = THEMES[theme]
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter, NullLocator
    fmt = FuncFormatter(lambda v, _: ("%g" % v) if 1e-3 <= v <= 1e4 else (r"$10^{%d}$" % round(math.log10(v))))
    for name, key, xs, sc, xl, xt in PANELS:
        h, e, tot = sweep(key, xs)
        fig, ax = new_fig()
        style(ax, t)
        ax.plot(xs * sc, h, color=t["hyst"], lw=1.6, label="magnetic")
        ax.plot(xs * sc, e, color=t["eddy"], lw=1.6, label="eddy current")
        ax.plot(xs * sc, tot, color=t["total"], lw=2.2, label="total")
        x0 = BASE[key] * sc
        ax.axvline(x0, color=t["muted"], lw=0.8, ls=(0, (3, 3)))
        ax.plot([x0], [1.0], "o", ms=5, color=t["total"])
        if name == "K":
            # below this K_eff the working frequency is above f_r: the dynamic term does not apply there
            Kx = math.pi * BASE["Js"] * BASE["f"] / (GAMMA * MU0) * sc
            ax.axvspan(xs[0] * sc, Kx, color=t["grid"], alpha=0.5, lw=0)
            ax.axvline(Kx, color=t["muted"], lw=0.8, ls=":")
            ax.text(Kx, 0.04, r"$f_\mathrm{r}$ ", transform=ax.get_xaxis_transform(), color=t["muted"],
                    fontsize=9, ha="right", va="bottom")
        if name == "f":
            xr = f_r(BASE) * sc
            ax.axvline(xr, color=t["muted"], lw=0.8, ls=":")
            ax.text(xr, 0.04, r" $f_\mathrm{r}$", transform=ax.get_xaxis_transform(), color=t["muted"],
                    fontsize=9, ha="left", va="bottom")
        if name == "d":
            xo = d_opt(BASE) * sc
            ax.axvline(xo, color=t["muted"], lw=0.8, ls=":")
            ax.text(xo, 0.04, r" $d_\mathrm{opt}$", transform=ax.get_xaxis_transform(), color=t["muted"],
                    fontsize=9, ha="left", va="bottom")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_locator(FixedLocator(xt))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: r"$10^{%d}$" % round(math.log10(v))) if name == "mur"
                                     else fmt)
        span = math.log10(max(h.max(), e.max(), tot.max()) / min(h.min(), e.min(), tot.min()))
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,) if span > 2.5 else (1.0, 2.0, 5.0)))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(fmt)
        ax.set_xlabel(xl, fontsize=10)
        ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
        lg = ax.legend(frameon=False, fontsize=8.5, loc="best")
        for tx in lg.get_texts():
            tx.set_color(t["fg"])
        save(fig, name, theme)


def homogeneity_panel(theme):
    t = THEMES[theme]
    fig, ax = new_fig()
    style(ax, t)
    hs = np.linspace(0.0, 1.0, 101)
    h, e, tot = sweep("h", hs)
    ax.plot(hs, h, color=t["hyst"], lw=1.6, label="magnetic")
    ax.plot(hs, e, color=t["eddy"], lw=1.6, label="eddy current")
    ax.plot(hs, tot, color=t["total"], lw=2.2, label="total")
    ax.plot([BASE["h"]], [1.0], "o", ms=5, color=t["total"])
    ax.set_xlabel(r"inhomogeneity $h$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    lg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    for tx in lg.get_texts():
        tx.set_color(t["fg"])
    save(fig, "h", theme)


def contact_panels(theme):
    t = THEMES[theme]
    # loss against the contact probability, one line per contact area fraction
    fig, ax = new_fig()
    style(ax, t)
    ps = np.linspace(0.0, 1.0, 101)
    for col, a in zip(t["ramp"], (0.0, 0.05, 0.1, 0.2)):
        _, _, tot = sweep("p_c", ps, a_c=a)
        ax.plot(ps, tot, color=col, lw=1.8, label=r"$a_\mathrm{c}$ = %.2f" % a)
    ax.plot([BASE["p_c"]], [1.0], "o", ms=5, color=t["total"])
    ax.set_xlabel(r"contact probability $p_\mathrm{c}$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 1)
    lg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    for tx in lg.get_texts():
        tx.set_color(t["fg"])
    save(fig, "pc", theme)
    # loss against the contact area fraction, one line per contact probability
    fig, ax = new_fig()
    style(ax, t)
    acs = np.linspace(0.0, 0.3, 101)
    for col, pc in zip(t["ramp"], (0.0, 0.2, 0.5, 1.0)):
        _, _, tot = sweep("a_c", acs, p_c=pc)
        ax.plot(acs, tot, color=col, lw=1.8, label=r"$p_\mathrm{c}$ = %.1f" % pc)
    ax.plot([BASE["a_c"]], [1.0], "o", ms=5, color=t["total"])
    ax.set_xlabel(r"contact area fraction $a_\mathrm{c}$", fontsize=10)
    ax.set_ylabel(r"$P\,/\,P_\mathrm{base}$", fontsize=10)
    ax.set_xlim(0, 0.3)
    lg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    for tx in lg.get_texts():
        tx.set_color(t["fg"])
    save(fig, "ac", theme)


if __name__ == "__main__":
    plt.rcParams["mathtext.fontset"] = "cm"
    for th in THEMES:
        power_panels(th)
        homogeneity_panel(th)
        contact_panels(th)
    ph, pe, fc = parts(BASE)
    print("base: hysteresis %.3g W/m3, eddy %.3g W/m3, contact factor %.3f, P_base %.3g W/m3 = %.3g W/cm3"
          % (ph, pe, fc, P0, P0 * 1e-6))
    print("mu_r %.3g; f_r %.0f kHz; share hysteresis %.2f, eddy %.2f; phi_s %.4f; d_opt %.2f um (thin-layer form %.2f um)" % (
        mu_r(BASE), f_r(BASE) / 1e3, ph / (ph + pe), pe / (ph + pe), phi_s(BASE), 1e6 * d_opt(BASE), 1e6 * d_opt_thin(BASE)))
