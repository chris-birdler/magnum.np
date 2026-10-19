# Audit: current values of the deck numbers

Deck built 11 Oct 2026 16:45 from results/deck_data.json (11 Oct 16:42, analysis file sobol_analysis_steady2.json, reruns used in the run list: 0). Current = audit recomputation with the main rule of PROTOKOLL 7.65.

- β all: deck 1.56 … 2.60, mean 2.11 (sd 0.25), > 2: 78 of 116 → current 1.59 … 2.60, mean 2.116 (sd 0.250), > 2: 77 of 116
- repeat scatter of β all (pooled, 20 repeats in 5 groups, dof 15): deck 0.14 → current 0.131; per group: C 0.152 (n 8), K1 0.115 (n 3), K2 0.119 (n 3), K3 0.085 (n 3), K4 0.116 (n 3)
- spread from the parameters √(sd² − sd_rep²): deck 0.21 → current 0.212
- loop-to-loop scatter of P (median sd of ln P over calm tails with ≥ 3 loops): deck 2.5 % → current 2.50 % (n 783 stages)

Main fit (audit recomputation, identical to results/sobol_analysis_steady2_mc2.json to 4 digits):

- β × lnQ (deck sign: A up = lnd down): +0.150 (+0.065 … +0.227), p 0.0002, class lever
- β × r_p (deck sign: A up = lnd down): +0.145 (+0.027 … +0.258), p 0.0135, class not determined
- β × phi (deck sign: A up = lnd down): -0.064 (-0.182 … +0.024), p 0.2033, class not determined
- β × lnd (deck sign: A up = lnd down): +0.572 (+0.479 … +0.663), p 0.0000, class lever
- P@50 × lnQ: -6.3 % (-22 … +12 %), class not determined
- P@50 × r_p: +3.7 % (-18 … +32 %), class not determined
- P@50 × phi: +6.6 % (-13 … +32 %), class not determined
- P@50 × lnd: -26.5 % (-40 … -12 %), class lever
- A effect on ln P per amplitude (A up): 9 mT -1.12, 20 mT -1.01, 35 mT -0.69, 50 mT -0.31, 70 mT +0.01, 100 mT +0.19, 150 mT +0.28; sign change at ≈ 70 mT (deck: ≈ 69 mT)

Mesh pair K3 (A = 20) against base (A = 10) (K3 also differs in L_eff, r_p, φ):

- dx 1.5: 9 mT × 0.82 (ln -0.20 ± 0.18, p 0.308, n 6 / 4); 50 mT × 1.59 (ln +0.47 ± 0.18, p 0.034, n 6 / 4); 100 mT × 1.91 (ln +0.65 ± 0.08, p 0.000, n 6 / 4); crossing ≈ 15 mT
- dx 3: 9 mT × 0.35 (ln -1.06 ± 0.30, p 0.029, n 6 / 4); 50 mT × 0.79 (ln -0.23 ± 0.15, p 0.152, n 6 / 4); 100 mT × 1.34 (ln +0.29 ± 0.09, p 0.010, n 6 / 4); crossing ≈ 68 mT
- Δβ (9/50/100 fit) K3 − base, dx 3: +0.539 ± 0.110 (n 6 / 4, Welch p 0.010)
- Δβ (9/50/100 fit) K3 − base, dx 1.5: +0.363 ± 0.066 (n 6 / 4, Welch p 0.002)
- K3 β 9 … 100 mT: dx 3 2.372, dx 1.5 2.302 (means of per-seed fits), paired +0.071 ± 0.029 (n 6, p 0.057)
