import torch
import os
from magnumnp.common import logging, DecoratedTensor

__all__ = ["State"]

class State(object):
    def __init__(self, mesh, t0=0, device=None):
        self.t = t0
        self.mesh = mesh
        if device == None:
            CUDA_DEVICE = os.environ.get('CUDA_DEVICE', '0')
            self._device = torch.device(f"cuda:{CUDA_DEVICE}" if torch.cuda.is_available() else "cpu")
        else:
            self._device = device
        self.material = {}
        logging.info_green("[State] running on device:%s" % self._device)
        logging.info_green("[Mesh] %dx%dx%d (size= %g x %g x %g)" % (mesh.n + mesh.dx))

    def _zeros(self, size, dtype=torch.float64, **kwargs):
        return torch.zeros(size, dtype=dtype, device=self._device, **kwargs)

    def _arange(self, start, end=None, step=1, dtype=torch.float64, **kwargs):
        if end == None:
           end = start
           start = 0
        return torch.arange(start, end, step, dtype=dtype, device=self._device, **kwargs)

    def _linspace(self, start, end, steps, dtype=torch.float64, **kwargs):
        return torch.linspace(start, end, steps, dtype=dtype, device=self._device, **kwargs)

    # TODO: make this _tensor
    def tensor(self, data, dtype=torch.float64):
        return torch.tensor(data, dtype=dtype, device=self._device)

    def Tensor(self, data, dtype=torch.float64, requires_grad = False):
        if isinstance(data, list) or isinstance(data, tuple):
            t = torch.tensor(data, dtype=dtype, device=self._device).as_subclass(DecoratedTensor)
            t.requires_grad = requires_grad
            return t
        elif isinstance(data, torch.Tensor):
            return data.as_subclass(DecoratedTensor)
        elif callable(data): # TODO: allow returning List instead of Tensors
            return data
        else:
            raise TypeError("h needs to be 'list', 'torch.Tensor', or 'function'!")

    def Constant(self, c):
        c = self.tensor(c)
        x = self._zeros(self.mesh.n + c.shape).as_subclass(DecoratedTensor)
        x[...] =  c
        return x

    def SpatialCoordinates(self):
        x = self._arange(self.mesh.n[0]) * self.mesh.dx[0] + self.mesh.dx[0]/2. + self.mesh.origin[0]
        y = self._arange(self.mesh.n[1]) * self.mesh.dx[1] + self.mesh.dx[1]/2. + self.mesh.origin[1]
        z = self._arange(self.mesh.n[2]) * self.mesh.dx[2] + self.mesh.dx[2]/2. + self.mesh.origin[2]

        XX, YY, ZZ = torch.meshgrid(x, y, z, indexing = "ij")
        return XX, YY, ZZ
