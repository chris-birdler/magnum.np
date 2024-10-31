# %% [markdown]
# This jupyter-notebook has be created with 'jupytext'

# %% [markdown]
# ## Install magnum.np and fetch reference data (required for Colab)

# %%
!pip install -q triton magnumnp
from os import path
if not path.isdir("ref"):
    !mkdir ref
    !wget -P ref https://gitlab.com/magnum.np/magnum.np/raw/main/demos/sp4/ref/m.dat &> /dev/null
