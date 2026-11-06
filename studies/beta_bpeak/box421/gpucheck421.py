# GPU load check before the set-up: fp32 matmul throughput and device copy bandwidth per GPU.
# Free V100: about 13 TFLOPS and about 800 GB/s.
import time
import torch

for d in range(torch.cuda.device_count()):
    dev = torch.device("cuda", d)
    a = torch.randn(8192, 8192, device=dev)
    b = torch.randn(8192, 8192, device=dev)
    for _ in range(3):
        a @ b
    torch.cuda.synchronize(dev)
    t = time.time()
    n = 20
    for _ in range(n):
        a @ b
    torch.cuda.synchronize(dev)
    tf = 2 * 8192 ** 3 * n / (time.time() - t) / 1e12
    x = torch.empty(256 * 2 ** 20, device=dev)  # 1 GiB
    y = torch.empty_like(x)
    for _ in range(3):
        y.copy_(x)
    torch.cuda.synchronize(dev)
    t = time.time()
    for _ in range(n):
        y.copy_(x)
    torch.cuda.synchronize(dev)
    bw = 2 * x.numel() * 4 * n / (time.time() - t) / 1e9
    print("gpu %d %s: %.1f TFLOPS fp32, %.0f GB/s copy" % (d, torch.cuda.get_device_name(d), tf, bw))
