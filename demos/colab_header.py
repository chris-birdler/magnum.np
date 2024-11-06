# %% [markdown]
# This jupyter-notebook has be created with 'jupytext'.
#
# Use this [Link](https://colab.research.google.com/github/magnumnp/magnumnp_demos/blob/main/DEMO.ipynb) to directly open the Notebook in Google Colab.

# %% [markdown]
# ## Install magnum.np and fetch reference data (required for Colab)

# ---
# jupyter:
#   metadata:
#     colab:
#       gpuType: T4
#     accelerator: GPU
# ---

# %%
!pip install -q triton magnumnp
from os import path
if not path.isdir("ref"):
    !mkdir ref
    !wget -P ref https://gitlab.com/magnum.np/magnum.np/raw/main/demos/DEMO/ref/m.dat &> /dev/null
