"""Build loss_equation.html: the showroom loss equation with the trend plots of plots.py embedded as data URIs.

Chris 2026-11-04: no verbal statements on the page. It shows only the equations, the meaning and the base value of
each variable, and the plots. All values come from the model in plots.py.
"""
import base64
import pathlib

import plots as m

HERE = pathlib.Path(__file__).resolve().parent
B = m.BASE


def img(name):
    out = []
    for th in ("light", "dark"):
        b = base64.b64encode((HERE / ("fig_%s_%s.png" % (name, th))).read_bytes()).decode()
        out.append('<img class="fig-%s" src="data:image/png;base64,%s" alt="">' % (th, b))
    return "\n".join(out)


def g(v, digits=3):
    """A number with at most `digits` significant digits, without trailing zeros."""
    return "%.*g" % (digits, v)


CARDS = [
    ("h", "Inhomogeneity", r"P_\mathrm{m} \propto 1 + \gamma_\mathrm{h}\, h"),
    ("K", "Effective anisotropy",
     r"P_\mathrm{m} \propto K_\mathrm{eff}\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad f_\mathrm{r} \propto K_\mathrm{eff}"),
    ("Js", "Saturation polarization",
     r"P_\mathrm{m} \propto J_\mathrm{s}^{-\beta}\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad f_\mathrm{r} \propto J_\mathrm{s}^{-1}"),
    ("B", "Peak flux density", r"P_\mathrm{m} \propto B_\mathrm{peak}^{\beta}, \quad P_\mathrm{e} \propto B_\mathrm{peak}^{2}"),
    ("f", "Frequency", r"P_\mathrm{m} \propto f\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad P_\mathrm{e} \propto f^{2}"),
    ("d", "Particle diameter", r"P_\mathrm{m} \propto 1 + \gamma_\mathrm{s}\,\varphi_\mathrm{s}(d), \quad P_\mathrm{e} \propto d^{2}"),
    ("rho", "Resistivity", r"P_\mathrm{e} \propto \rho^{-1}"),
    ("pc", "Contact probability", r"P \propto 1 + p_\mathrm{c}\,(\kappa_0 + \kappa_\mathrm{a} a_\mathrm{c})"),
    ("ac", "Contact area", r"P \propto 1 + p_\mathrm{c}\,(\kappa_0 + \kappa_\mathrm{a} a_\mathrm{c})"),
]

# symbol, meaning, unit, base value
INPUTS = [
    (r"f", "frequency", "Hz", "%s MHz" % g(B["f"] / 1e6)),
    (r"B_\mathrm{peak}", "peak flux density of the core", "T", "%s mT" % g(B["B"] * 1e3)),
    (r"J_\mathrm{s}", "saturation polarization of the particle material", "T", "%s T" % g(B["Js"])),
    (r"K_\mathrm{eff}", "effective anisotropy: crystal anisotropy and magnetostriction", "J/m³",
     "%s kJ/m³" % g(B["K"] / 1e3)),
    (r"h", "inhomogeneity: anisotropy fluctuation between regions or particles, relative to \\(K_\\mathrm{eff}\\)",
     "–", g(B["h"])),
    (r"d", "particle diameter", "m", "%s µm" % g(B["d"] * 1e6)),
    (r"\delta_\mathrm{s}", "thickness of the disturbed surface layer of a particle", "m", "%s µm" % g(B["ds"] * 1e6)),
    (r"\rho", "resistivity of the particle material", "Ω·m", "%s µΩ·m" % g(B["rho"] * 1e6)),
    (r"p_\mathrm{c}", "contact probability: 0 all particles insulated, 1 all particles in contact", "–", g(B["p_c"])),
    (r"a_\mathrm{c}", "contact area per particle surface", "–", g(B["a_c"])),
]

# symbol, meaning, value
CONSTS = [
    (r"\beta", "flux density exponent of the magnetic part", g(B["beta"])),
    (r"\alpha", "frequency exponent of the dynamic part", g(B["alpha"])),
    (r"c_\mathrm{m}", "shape factor of the magnetic loss per cycle", g(B["c_h"])),
    (r"\gamma_\mathrm{h}", "loss increase per unit of inhomogeneity", g(B["gh"])),
    (r"\gamma_\mathrm{s}", "extra loss of the surface layer, relative to the particle core", g(B["gs"])),
    (r"\kappa_0", "loss increase per contact probability", g(B["k0"])),
    (r"\kappa_\mathrm{a}", "extra slope per unit of \\(a_\\mathrm{c}\\)", g(B["ka"])),
    (r"\gamma", "gyromagnetic ratio", "%s·10¹¹ s⁻¹ T⁻¹" % g(m.GAMMA / 1e11)),
    (r"\mu_0", "magnetic constant", "4π·10⁻⁷ T·m/A"),
]

PM, PE, FC = m.parts(B)
# symbol, meaning, unit, value at the base point
DERIVED = [
    (r"P", "loss per volume of particle material", "W/m³", "–"),
    (r"P_\mathrm{m}", "magnetic part of \\(P\\), without the contact factor", "W/m³", "%s W/cm³" % g(PM * 1e-6, 2)),
    (r"P_\mathrm{e}", "eddy current part of \\(P\\), without the contact factor", "W/m³", "%s W/cm³" % g(PE * 1e-6, 2)),
    (r"P_\mathrm{base}", "\\(P\\) at the base point", "W/m³", "%s W/cm³" % g(m.P0 * 1e-6, 2)),
    (r"\mu_\mathrm{r}", "relative permeability of the particle material", "–", g(m.mu_r(B))),
    (r"f_\mathrm{r}", "resonance frequency: there the dynamic part is as large as the static part", "Hz",
     "%s MHz" % g(m.f_r(B) / 1e6)),
    (r"\varphi_\mathrm{s}", "volume fraction of the surface layer in a particle", "–", g(m.phi_s(B), 2)),
    (r"S/V", "surface-to-volume ratio of a particle", "1/m", "%s·10⁶ 1/m" % g(6.0 / B["d"] / 1e6)),
    (r"d_\mathrm{opt}", "particle diameter with the smallest loss", "m", "%s µm" % g(m.d_opt(B) * 1e6, 2)),
    (r"K_1", "crystal anisotropy of a grain", "J/m³", "–"),
    (r"D", "grain size", "m", "–"),
    (r"A", "exchange constant", "J/m", "–"),
    (r"L_0", "exchange length of the grain anisotropy", "m", "–"),
]


def rows(items):
    return "\n".join("<tr>%s</tr>" % "".join(["<td>\\(%s\\)</td>" % it[0]] + ["<td>%s</td>" % x for x in it[1:]])
                     for it in items)


cards = "\n".join(
    '<figure class="card" id="%s"><figcaption><h3>%s</h3><div class="law">\\[%s\\]</div></figcaption>'
    '<div class="plot">%s</div></figure>' % (n, h, law, img(n)) for n, h, law in CARDS)

page = (HERE / "template.html").read_text()
page = (page.replace("{{CARDS}}", cards).replace("{{INPUTS}}", rows(INPUTS)).replace("{{CONSTS}}", rows(CONSTS))
        .replace("{{DERIVED}}", rows(DERIVED)))
assert "{{" not in page
(HERE / "loss_equation.html").write_text(page)
print("written %d kB" % (len(page) // 1024))
