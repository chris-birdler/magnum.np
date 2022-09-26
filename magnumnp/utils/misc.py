import torch

__all__ = ["add_noise", "nsk"]

def add_noise(x, dev = 1.0, mean = 0.0):
   if torch.is_tensor(x):
        x += torch.empty_like(x).normal_(mean = mean, std = dev)
        x.normalize()

def nsk(x, axis = 2):
    m = self.mean(axis=axis)
    dxm = torch.stack(torch.gradient(m, spacing = dx[0], dim = 0), dim = -1).squeeze(-1)
    dym = torch.stack(torch.gradient(m, spacing = dx[1], dim = 1), dim = -1).squeeze(-1)
    nsk = 1./(4. * pi) * (m * torch.cross(dxm, dym)).sum()* dx[0] * dx[1]
    return float(nsk.detach().cpu().numpy())

