# Showroom loss equation

A one-page equation for the trends of the power loss of soft magnetic powder cores (request of Chris, 2026-11-04;
artifact "Powder Core Loss Equation", https://claude.ai/artifact/1RkYUK4t2M5i1FowmbXS9H, version 11). The deck with the
same equation, base point and plots is in the repository slide_decks, folder showroom_loss_equation.

Since version 11 (Chris 2026-11-04): no prose, only the equations, the meaning and base value of each variable and the
plots. Base point for the MHz comparison: f 30 MHz (plots 1 to 50 MHz), B_peak 10 mT, J_s 1.5 T, K_eff 20 kJ/m3
(mu_r 30.8, f_r 939 MHz), d 3 um, rho 0.43 uOhm m (iron-nickel-like), c_m 0.048 and gamma_s 0.26 (4.4 W/cm3 at the base
point, d_opt 3.3 um at 1 MHz). Version 10 (100 kHz, nanocrystalline base point, text cards, limits, sources) is in the
git history (commit 702cf15).

- plots.py: the model, the base point and the trend plots (fig_<name>_light.png and fig_<name>_dark.png).
- build_page.py, template.html: build loss_equation.html (the plots embedded as data URIs).

    python plots.py && python build_page.py
