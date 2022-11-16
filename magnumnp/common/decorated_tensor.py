import torch

__all__ = ["DecoratedTensor"]

class DecoratedTensor(torch.Tensor):
    @staticmethod
    def __new__(cls, x, *args, **kwargs): # TODO: is this needed?
        return super().__new__(cls, x, *args, **kwargs)

    def avg(self, dim=(0,1,2)):
        if self.dim() <= 1: # e.g. [0,0,1]
            return self
        elif self.dim() == 2: # state.m[domain]
            return self.mean(dim=0)
        else:
            return self.mean(dim=dim)

    def average(self, dim=(0,1,2)):
        return self.avg(dim)

    def normalize(self):
        self /= torch.linalg.norm(self, dim = -1, keepdim = True)
        self[...] = torch.nan_to_num(self, posinf=0, neginf=0)
        return self

    def pad(self, dim, n):
        shape = list(self.shape)
        shape[dim] = abs(n)
        zeros = torch.zeros(shape)

        if n > 0:
            return torch.concat([self, zeros], dim=dim)
        elif n < 0:
            return torch.concat([zeros, self], dim=dim)
        else:
            return self

    def __call__(self, t):
        return self

    def __getitem__(self, key):
        if self.dim() == 0 or (self.dim() == 1 and self.shape[0] == 1):
            return self
        else:
            super().__getitem__(key)
