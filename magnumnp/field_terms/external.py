import torch
from scipy import constants
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

__all__ = ["ExternalField"]

class ExternalField(object):
    def __init__(self, mesh, material, h):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._h = torch.zeros(mesh.n + (3,), dtype=torch.float64, device = cuda)
        self._h[:,:,:,:] = torch.DoubleTensor(h)

    def h(self, t, m):
        return self._h

    def E(self, t, m):
        return - constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._Ms * m * self.h(t, m))
