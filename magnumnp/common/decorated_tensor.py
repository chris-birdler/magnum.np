import torch

__all__ = ["DecoratedTensor"]

class DecoratedTensor(torch.Tensor):
    @staticmethod 
    def __new__(cls, x, *args, **kwargs): 
        return super().__new__(cls, x, *args, **kwargs) 
      
    def avg(self, dim=(0,1,2)):
        if self.dim() == 1: # e.g. [0,0,1]
            return self
        elif self.dim() == 2: # state.m[domain]
            return self.mean(dim=0) 
        else:
            return self.mean(dim=dim) 

    def normalize(self):
        self /= torch.linalg.norm(self, dim = 3, keepdim = True)
        self[...] = torch.nan_to_num(self, posinf=0, neginf=0)

    def __call__(self, t):
        return self
