"""Build loss_equation.html: the showroom loss equation with the trend plots of plots.py embedded as data URIs."""
import base64
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def img(name):
    out = []
    for th in ("light", "dark"):
        b = base64.b64encode((HERE / ("fig_%s_%s.png" % (name, th))).read_bytes()).decode()
        out.append('<img class="fig-%s" src="data:image/png;base64,%s" alt="">' % (th, b))
    return "\n".join(out)


CARDS = [
    ("h", "Homogeneity", r"P_\mathrm{m} \propto 1 + \gamma_\mathrm{h}\, h",
     "\\(h\\) is the relative fluctuation of the anisotropy between regions or particles (stress, composition), "
     "on scales larger than the exchange length. A perfectly homogeneous material has "
     "\\(h\\) of 0 and the smallest loss. The magnetic part rises linearly with \\(h\\)."),
    ("K", "Effective anisotropy", r"P_\mathrm{m} \propto K_\mathrm{eff}\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad f_\mathrm{r} \propto K_\mathrm{eff}",
     "\\(K_\\mathrm{eff}\\) contains all anisotropy: the crystal anisotropy (for nanocrystalline materials after "
     "the Herzer averaging) and the magnetoelastic part. The static part is proportional to it. A higher "
     "\\(K_\\mathrm{eff}\\) also means a lower permeability and a higher resonance frequency, so the dynamic part "
     "rises more slowly. In the grey range the working frequency is above \\(f_\\mathrm{r}\\); the dynamic term "
     "does not apply there."),
    ("Js", "Saturation polarization", r"P_\mathrm{m} \propto J_\mathrm{s}^{-\beta}\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad f_\mathrm{r} \propto J_\mathrm{s}^{-1}",
     "At the same flux density, a material with a higher \\(J_\\mathrm{s}\\) is driven less far towards "
     "saturation, so the magnetic part falls. Its permeability is higher and its resonance lower, which slows the "
     "fall a little. The eddy current part does not change."),
    ("B", "Peak flux density", r"P_\mathrm{m} \propto B_\mathrm{peak}^{\beta}, \quad P_\mathrm{e} \propto B_\mathrm{peak}^{2}",
     "Both parts rise with the drive amplitude. With the base exponent of 2, the two lines are parallel."),
    ("f", "Frequency", r"P_\mathrm{m} \propto f\,\big[1 + (f/f_\mathrm{r})^{\alpha-1}\big], \quad P_\mathrm{e} \propto f^{2}",
     "Far below \\(f_\\mathrm{r}\\) (dotted line) the magnetic part is the static hysteresis, a fixed loss per "
     "cycle. Towards \\(f_\\mathrm{r}\\) the dynamic part adds, and the slope rises from 1 towards "
     "\\(\\alpha\\). The eddy current part rises with the square of the frequency. The plot stops at about twice "
     "\\(f_\\mathrm{r}\\), because the dynamic term applies only below and near \\(f_\\mathrm{r}\\)."),
    ("d", "Particle diameter", r"P_\mathrm{m} \propto 1 + \gamma_\mathrm{s}\,\varphi_\mathrm{s}(d), \quad P_\mathrm{e} \propto d^{2}",
     "Small particles increase the magnetic part, because the disturbed surface layer becomes a larger part of "
     "each particle; in very small particles the whole volume is surface layer and the loss stops rising. Large "
     "particles increase the eddy current part. The total loss has its minimum at the optimum diameter (dotted "
     "line): small particles are not always better."),
    ("rho", "Resistivity", r"P_\mathrm{e} \propto \rho^{-1}",
     "A higher resistivity of the particle material decreases the eddy current part."),
    ("pc", "Contact probability", r"P \propto 1 + p_\mathrm{c}\,(\kappa_0 + \kappa_a a_\mathrm{c})",
     "The loss rises linearly from a fully insulated powder (0) to a powder in which all particles touch (1). "
     "Each line is one contact area fraction: a larger contact area makes the line steeper."),
    ("ac", "Contact area", r"P \propto 1 + p_\mathrm{c}\,(\kappa_0 + \kappa_a a_\mathrm{c})",
     "At a given contact probability, the loss rises linearly with the contact area. Without contacts (bottom "
     "line) the contact area has no effect."),
]

INPUTS = [
    (r"f", "frequency", "Hz", "100 kHz", "showroom choice"),
    (r"B_\mathrm{peak}", "peak flux density of the core (amplitude)", "T", "100 mT", "showroom choice"),
    (r"J_\mathrm{s}", "saturation polarization, \\(\\mu_0 M_\\mathrm{s}\\)", "T", "1.25 T",
     "typical of nanocrystalline FeSiBNbCu"),
    (r"K_\mathrm{eff}", "effective anisotropy: crystal anisotropy (after the Herzer averaging for nanocrystalline "
     "materials) and magnetostriction", "J/m³", "4 J/m³",
     "Herzer average for stress-free FeSiBNbCu grains of 10 nm (box below): a lower bound; stress in a powder can "
     "raise it to about 100 J/m³ (estimate of our record)"),
    (r"h", "inhomogeneity: relative fluctuation of the anisotropy, \\(\\Delta K / K_\\mathrm{eff}\\), on scales "
     "larger than the exchange length (between regions or particles: stress, composition); 0 is perfectly "
     "homogeneous", "–", "0.2", "showroom choice"),
    (r"d", "particle diameter", "m", "10 µm", "showroom choice"),
    (r"\delta_\mathrm{s}", "thickness of the disturbed surface layer of a particle (stress, oxide, defects)", "m",
     "0.1 µm", "illustrative"),
    (r"\rho", "resistivity of the particle material", "Ω·m", "1.2 µΩ·m", "measured 1.17 µΩ·m [3]"),
    (r"p_\mathrm{c}", "contact probability: 0 all particles insulated, 1 all particles in contact", "–", "0.2",
     "showroom choice"),
    (r"a_\mathrm{c}", "contact area per particle surface", "–", "0.05", "showroom choice"),
]

CONSTS = [
    (r"\beta", "flux density exponent of the magnetic part", "2",
     "Steinmetz β usually 2 to 3 [8]; a powder core datasheet 1.6 [9]; Rayleigh range 3. Our micromagnetic "
     "model: 1.9 to 2.3 for the total rate-dependent loss at 30 MHz (4 design points, \\(K_\\mathrm{eff}\\) 1 to "
     "6 kJ/m³, coarse mesh)"),
    (r"c_\mathrm{m}", "shape factor of the magnetic loss per cycle, order 1", "1", "illustrative"),
    (r"\alpha", "frequency exponent of the dynamic part", "1.5",
     "shape borrowed from the excess loss of wall motion, 1.5 [10, 12]; a linear damping relaxation gives at most "
     "2. For "
     "comparison, the Steinmetz α of the total loss: usually 1 to 2 [8], a powder core datasheet 1.0 below "
     "100 kHz and about 1.4 above [9]"),
    (r"\gamma", "gyromagnetic ratio", "1.76·10¹¹ s⁻¹T⁻¹", "physical constant"),
    (r"\mu_0", "magnetic constant", "4π·10⁻⁷ T·m/A", "physical constant"),
    (r"\gamma_\mathrm{h}", "loss increase per unit of inhomogeneity", "2", "illustrative, set by the homogeneity rule"),
    (r"\gamma_\mathrm{s}", "additional loss of the surface layer relative to the particle core", "1.2",
     "illustrative; chosen so that the optimum particle diameter is 3 to 4 µm (Chris, 2026-11-04)"),
    (r"\kappa_0", "loss increase per contact probability", "1", "illustrative, set by the contact rule"),
    (r"\kappa_a", "additional slope per contact area fraction", "20", "illustrative, set by the contact rule"),
]


def rows(items, n):
    out = []
    for it in items:
        cells = ["<td>\\(%s\\)</td>" % it[0]] + ["<td>%s</td>" % x for x in it[1:n]]
        out.append("<tr>%s</tr>" % "".join(cells))
    return "\n".join(out)


cards = "\n".join(
    '<figure class="card" id="%s"><figcaption><h3>%s</h3><div class="law">\\[%s\\]</div><p>%s</p></figcaption>'
    '<div class="plot">%s</div></figure>' % (n, h, law, txt, img(n)) for n, h, law, txt in CARDS)

page = (HERE / "template.html").read_text()
page = page.replace("{{CARDS}}", cards).replace("{{INPUTS}}", rows(INPUTS, 5)).replace("{{CONSTS}}", rows(CONSTS, 4))
(HERE / "loss_equation.html").write_text(page)
print("written %d kB" % (len(page) // 1024))
