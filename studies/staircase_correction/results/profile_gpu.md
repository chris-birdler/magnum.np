# GPU profile of the staircase correction (PROTOKOLL 3.35)

Tesla V100-SXM2-16GB (vast.ai 57100794), profile_staircase.py (commit e48d6d8), single precision, CUDA_DEVICE 0, 20 repetitions after a warm-up. Times in ms per field evaluation. Base: Newell stray field + exchange + anisotropy. GPU memory: peak of the set-up of the correction, allocated / reserved (GB). Test geometry: T3Pc_s941 with its m; the other geometries: random m.

| case | N | d/l_ex | q | pairs (M) | cut cells | base | correction | pair sum | centroid term | factor of the field terms | set-up (s) | GPU memory (GB) | host peak (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| test_q2 | 48 | 92 | 2 | 2.7 | 8736 | 1.8 | 7.2 | 2.2 | 6.4 | 5.0 | 1.8 | 0.64 / 1.13 | 1.4 |
| test_q4 | 48 | 92 | 4 | 4.3 | 13760 | 1.7 | 7.0 | 3.4 | 6.3 | 5.1 | 2.0 | 0.94 / 1.64 | 1.7 |
| test_q8 | 48 | 92 | 8 | 4.9 | 15584 | 1.8 | 7.2 | 3.9 | 6.4 | 5.1 | 24.1 | 12.01 / 13.81 | 5.5 |
| d212_q2 | 104 | 211 | 2 | 15.3 | 46592 | 1.7 | 18.3 | 13.7 | 6.5 | 12.0 | 1.8 | 3.61 / 3.88 | 1.8 |
| d212_q4 | 104 | 211 | 4 | 23.1 | 70016 | 1.8 | 25.4 | 20.9 | 6.2 | 15.4 | 2.4 | 5.11 / 5.23 | 6.7 |
| d300_q2 | 148 | 301 | 2 | 31.1 | 93888 | 4.7 | 51.7 | 33.5 | 17.8 | 12.1 | 1.9 | 7.07 / 7.43 | 3.2 |
| d300_q4 | 148 | 301 | 4 | 47.4 | 142848 | 4.7 | 67.9 | 49.7 | 17.9 | 15.6 | 2.9 | 10.53 / 11.08 | 17.1 |
| d212_q8 | 104 | 211 | 8 | – | – | – | – | – | – | – | – | CUDA out of memory in the set-up (einsum, 4.02 GiB more at 15.66 GiB in use) | – |
| d300_q8 | 148 | 301 | 8 | – | – | – | – | – | – | – | – | – | process killed in the build of the fine geometry (1184³ cells) |
