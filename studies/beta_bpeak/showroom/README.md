# Showroom loss equation

A one-page equation for the trends of the power loss of soft magnetic powder cores (request of Chris, 2026-11-04;
artifact "Powder Core Loss Equation", version 10; reviewed twice by an independent agent). It shows directions and
strengths of trends, not loss values. The deck with the same equation is in the repository slide_decks,
folder showroom_loss_equation.

- plots.py: the model, the base point and the trend plots (fig_<name>_light.png and fig_<name>_dark.png).
- build_page.py, template.html: build loss_equation.html (the plots embedded as data URIs).

    python plots.py && python build_page.py
