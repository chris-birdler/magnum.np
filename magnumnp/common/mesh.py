import torch
import os

__all__ = ["Mesh"]

class Mesh(object):
    def __init__(self, n, dx, device=None):
        self.n = n
        self.dx = dx
        self.cell_volume = dx[0] * dx[1] * dx[2]
        if device == None:
            CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
            self.device = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
        else:
            self.device = device

    def function(self, dim):
        return torch.zeros(dim, dtype=torch.float64, device = self.device)

    def SpatialCoordinates(self):
        return x, y, z

    def __str__(self):
        return "%dx%dx%d_%gx%gx%g" % (self.n + self.dx)
