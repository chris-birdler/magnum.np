import torch
from scipy import constants
import os
CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
cuda = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")

__all__ = ["AnisotropyField"]

class AnisotropyField(object):
    def __init__(self, mesh, material):
        self._mesh = mesh
        self._Ms = material["Ms"]
        self._K = material["K"]
        self._K_axis = torch.zeros(mesh.n + (3,), dtype=torch.float64, device = cuda)
        self._K_axis[:,:,:,:] = torch.DoubleTensor(material["K_axis"])

        # initialize scratch space
        self._h = torch.zeros(mesh.n + (3,), dtype=torch.float64, device = cuda)

    def h(self, t, m):
        self._h[:,:,:,:] = 2. * self._K * self._K_axis / (constants.mu_0 * self._Ms) * torch.sum(self._K_axis * m, dim=3).unsqueeze(-1)
        self._h = torch.nan_to_num(self._h, posinf=0, neginf=0)
        return self._h

    def E(self, t, m):
        return - 0.5 * constants.mu_0 * self._mesh.cell_volume \
               * torch.sum(self._Ms * m * self.h(t, m))
