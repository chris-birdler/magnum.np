# AC demagnetization after the virgin relaxation: Bloch points and remanence (eval_predemag.py)

| run | variant | h0 | half cycles | BP virgin | BP max in the demag. | BP after relax 2 | remanence M/M_s | BP max at 9 mT | relax 2 iterations | demag. time (h) | ln w 9 mT − corner run |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K2_0_nx | 150 mT | 0.01547 | 19 | 9 | 10 | 1 | 0.0049 | 1 | 1500 | 1.69 | -0.490 |
| K2_1_nx | 150 mT | 0.01552 | 19 | 5 | 7 | 1 | 0.0011 | 1 | 13500 | 1.87 | +0.387 |
| K4_1_nx | 150 mT | 0.02004 | 21 | 2 | 2 | 0 | 0.0003 | 0 | 20000 | 2.09 | +0.253 |
| K3_1_nx | 150 mT | 0.00943 | 17 | 2 | 2 | 2 | 0.0022 | 2 | 7000 | 0.61 | -0.719 |

**Variant pd150 (150 mT): not clean.** K2_0_nx: 1 Bloch points after relax 2; K2_0_nx: Bloch points at 9 mT; K2_1_nx: 1 Bloch points after relax 2; K2_1_nx: Bloch points at 9 mT; K3_1_nx: 2 Bloch points after relax 2; K3_1_nx: Bloch points at 9 mT

| K2_0_nx | 2 × 150 mT | 0.03093 | 23 | 9 | 13 | 0 | 0.0058 | 0 | 10500 | 2.28 | -0.077 |
| K2_1_nx | 2 × 150 mT | 0.03103 | 23 | 5 | 10 | 0 | 0.0022 | 0 | 2500 | 1.68 | +0.067 |
| K4_1_nx | 2 × 150 mT | 0.04009 | 24 | 2 | 2 | 0 | 0.0004 | 0 | 7500 | 1.14 | +1.111 |
| K3_1_nx | 2 × 150 mT | 0.01886 | 20 | 2 | 2 | 2 | 0.0024 | 2 | 4500 | 0.43 | -0.521 |

**Variant pd300 (2 × 150 mT): not clean.** K3_1_nx: 2 Bloch points after relax 2; K3_1_nx: Bloch points at 9 mT

| K2_0_nx | saturation | 0.38 | 34 | 9 | 150 | 2 | 0.0058 | 2 | 10000 | 3.76 | +0.064 |
| K2_1_nx | saturation | 0.38 | 34 | 5 | 117 | 1 | 0.0012 | 1 | 10500 | 2.12 | +0.392 |
| K4_1_nx | saturation | 0.45 | 35 | 2 | 94 | 0 | 0.0006 | 0 | 4000 | 2.05 | -0.247 |
| K3_1_nx | saturation | 0.31 | 33 | 2 | 44 | 0 | 0.0009 | 0 | 1500 | 0.85 | -0.879 |

**Variant pdsat (saturation): not clean.** K2_0_nx: 2 Bloch points after relax 2; K2_0_nx: Bloch points at 9 mT; K2_1_nx: 1 Bloch points after relax 2; K2_1_nx: Bloch points at 9 mT

