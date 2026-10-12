"""Figure: beta (9 ... 150 mT) of the main-fit mean loss curve on a grid of A (d/l_ex at d = 1 µm) and K_eff (L_eff),
r_p and phi averaged over the design points; and P at 10 / 100 mT. dx 3 model (audit 2026-10-12)."""
import math

import numpy as np

import audit_lib as al
import design_data as dd
import sobol_audit as sa


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    recs = dd.load("steady2_mc2")
    X = dd.coded(recs)
    T = dd.terms(X)
    W = np.exp(np.array([r["lnw"] for r in recs]))
    vw, _, _ = sa.scatter_weights(recs, X)
    coefs = np.array([sa.irls_gamma(T, W[:, k], vw) for k in range(7)])
    g = np.linspace(-1, 1, 41)
    B = np.zeros((41, 41))
    P10 = np.zeros((41, 41))
    P100 = np.zeros((41, 41))
    pf = recs[0]["pf"]
    for i, xd in enumerate(g):
        for j, xq in enumerate(g):
            Xg = X.copy()
            Xg[:, 3], Xg[:, 0] = xd, xq
            lm = dd.terms(Xg) @ coefs.T                       # ln E[w] per design point and amplitude
            lmu = np.log(np.exp(lm).mean(0))                  # mean over r_p, phi of the design
            B[j, i] = lmu @ dd.OUT["beta_all"]
            P10[j, i] = math.exp(lmu @ dd.OUT["lnw@10mT"]) * pf
            P100[j, i] = math.exp(lmu @ dd.OUT["lnw@100mT"]) * pf
    lnd = dd.RANGES["lnd"][0] + (g + 1) / 2 * (dd.RANGES["lnd"][1] - dd.RANGES["lnd"][0])
    A = 895246.5549 * (1e-6 / np.exp(lnd)) ** 2 * 1e12               # pJ/m at d = 1 µm
    lnq = dd.RANGES["lnQ"][0] + (g + 1) / 2 * (dd.RANGES["lnQ"][1] - dd.RANGES["lnQ"][0])
    K = 895.2465549 * np.exp(lnq)                             # kJ/m^3
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, Z, lab, cm in ((axs[0], B, "β (9 … 150 mT) of the mean loss curve", "RdBu_r"),
                           (axs[1], P10, "P at 10 mT (W/cm³)", "viridis"), (axs[2], P100, "P at 100 mT (W/cm³)", "viridis")):
        cs = ax.contourf(A, K, Z, levels=14, cmap=cm)
        ax.contour(A, K, Z, levels=cs.levels, colors="k", linewidths=0.3)
        fig.colorbar(cs, ax=ax, label=lab)
        ax.set_xlabel("exchange stiffness A at d = 1 µm (pJ/m)\n(= particle size 1.0 → 0.71 µm at A = 10 pJ/m)")
        ax.set_ylabel("Herzer anisotropy K_eff (kJ/m³)")
        ax.set_yscale("log")
        for r in recs:
            ax.plot(r["A_pJm"], r["K_eff"] / 1e3, ".", color="0.3", ms=3)
    axs[0].set_title("dx 3 design model (Gamma GLM, 15 terms)\nr_p and φ averaged over the design points")
    fig.tight_layout()
    fig.savefig(al.HERE / "figures" / "map_beta_A_Keff.png", dpi=150)
    print("β corners: A20/K1.0 %.2f, A10/K1.0 %.2f, A20/K6.2 %.2f, A10/K6.2 %.2f" % (B[0, 0], B[0, -1], B[-1, 0], B[-1, -1]))
    print("A range", A.min(), A.max(), "K range", K.min(), K.max())


if __name__ == "__main__":
    main()
