import torch
import os
from magnumnp.common import logging

__all__ = ["State"]

class State(object):
    def __init__(self, mesh, material, t0=0, device=None):
        self.t = t0
        self._mesh = mesh
        self._material = material
        if device == None:
            CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
            self._device = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
        else:
            self._device = device
        logging.info_green("[State] running on device:%s" % self._device)

    def zeros(self, size, dtype=torch.float64, **kwargs):
        return torch.zeros(size, dtype=dtype, device=self._device, **kwargs)

    def arange(self, start, end=None, step=1, dtype=torch.float64, **kwargs):
        if end == None:
           end = start
           start = 0
        return torch.arange(start, end, step, dtype=dtype, device=self._device, **kwargs)

    def linspace(self, start, end, steps, dtype=torch.float64, **kwargs):
        return torch.linspace(start, end, steps, dtype=dtype, device=self._device, **kwargs)

    def DoubleTensor(self, data):
        return torch.DoubleTensor(data, device=self._device)

    def SpatialCoordinates(self):
        x = self.arange(self._mesh.dx[0]/2., (0.1 + self._mesh.n[0]) * self._mesh.dx[0], self._mesh.dx[0])
        y = self.arange(self._mesh.dx[1]/2., (0.1 + self._mesh.n[1]) * self._mesh.dx[1], self._mesh.dx[1])
        z = self.arange(self._mesh.dx[2]/2., (0.1 + self._mesh.n[2]) * self._mesh.dx[2], self._mesh.dx[2])

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return XX, YY, ZZ
