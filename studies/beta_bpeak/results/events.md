# Events in the measured cycles: physics or transient? (exploratory)

n = 115 α 0.1 runs × 7 amplitudes; event cycle: w_loop ≤ 0 or |w_dis / w_loop − 1| > 0.20.

## T1. Event rate against the cycle number (share of the run-amplitudes)

| amplitudes | n | cycle 1 | cycle 2 | cycle 3 | cycle 4 | cycle 5 | cycle 6 |
|---|---|---|---|---|---|---|---|
| 9 mT | 115 | 63.5 % | 35.7 % | 12.2 % | 3.5 % | 2.6 % | 0.0 % |
| 20 … 35 mT | 230 | 30.9 % | 13.0 % | 8.7 % | 6.1 % | 2.6 % | 1.7 % |
| 50 … 150 mT | 460 | 7.2 % | 3.5 % | 2.6 % | 1.5 % | 0.4 % | 0.0 % |

## T2. Energy sign of the event cycles (cycles 2 … 6) and ΣΔE

| amplitudes | event cycles | release (w_dis > w_loop) | uptake (w_dis < w_loop) | ΣΔE over cycles 2 … 6, median / loop loss (run-amplitudes with an event) | the same without an event |
|---|---|---|---|---|---|
| 9 mT | 62 | 60 | 2 | -0.80 (n 45) | -0.055 (n 70) |
| 20 … 35 mT | 74 | 74 | 0 | -2.43 (n 43) | -0.007 (n 187) |
| 50 … 150 mT | 37 | 37 | 0 | -1.10 (n 27) | -0.001 (n 433) |

## T3. Before and after an event (event in cycle 2 … 5; one cycle before and after)

| amplitudes | events | loop loss after / before: median | share lower after | n60 after − before: median | share with fewer n60 after | share with more n60 after |
|---|---|---|---|---|---|---|
| 9 mT | 60 | 0.82 | 67 % | -1 | 63 % | 17 % |
| 20 … 35 mT | 65 | 0.90 | 63 % | -2 | 65 % | 15 % |
| 50 … 150 mT | 37 | 0.76 | 78 % | -6 | 73 % | 16 % |

Reference, run-amplitudes without an event in cycles 1 … 6 (same cycle pairs): n60 after − before median +0, fewer 30 %, more 28 % (n 2344).

## T4. Closure of the calm cycles (cycles 2 … 6, no event)

| amplitudes | calm cycles | closure median | 90 % quantile | share > 1e-4 |
|---|---|---|---|---|
| 9 mT | 513 | 7.3e-05 | 2.0e-03 | 46 % |
| 20 … 35 mT | 1076 | 2.4e-04 | 2.1e-03 | 60 % |
| 50 … 150 mT | 2263 | 5.7e-04 | 2.1e-03 | 83 % |

## T5. Late events: closure of the 2 calm cycles before, against calm run-amplitudes

Late events: 12. Closure of the 2 calm cycles before: median 1.9e-03 (90 % quantile 2.8e-03). Calm run-amplitudes, cycles 4 … 6: median 3.2e-04 (90 % quantile 1.9e-03), n 2070.

Mann-Whitney (mean closure of the 2 cycles before each late event, n 12, against the calm cycles): p = 9.5e-06.

| run | B (mT) | cycle | closure 2 cycles before | 1 cycle before | balance of the event |
|---|---|---|---|---|---|
| S035 | 9 | 4 | 3.2e-03 | 1.6e-04 | 1.25 |
| K4_1 | 20 | 6 | 2.1e-03 | 3.2e-04 | 1.37 |
| S091 | 20 | 4 | 2.7e-03 | 7.8e-04 | 2.37 |
| C2 | 35 | 4 | 2.8e-03 | 2.1e-03 | 2.20 |
| S026 | 35 | 6 | 4.9e-03 | 2.1e-03 | 1.41 |
| S028 | 35 | 4 | 2.7e-03 | 2.3e-03 | 2.83 |
| S047 | 35 | 4 | 2.0e-03 | 1.2e-03 | 1.42 |
| S067 | 35 | 5 | 5.9e-04 | 1.1e-03 | 2.41 |
| S045 | 50 | 4 | 8.0e-04 | 1.6e-03 | 1.25 |
| S001 | 70 | 4 | 1.1e-03 | 2.0e-03 | 1.87 |
| K1_0 | 100 | 4 | 5.5e-04 | 1.3e-03 | 1.37 |
| S070 | 100 | 4 | 2.8e-03 | 1.7e-03 | 5.18 |

## T6. Mesh pairs: event rate in the cycles 1 … 4 at dx 3 and dx 1.5 (K3 and base, 9 / 50 / 100 mT)

| mesh | jobs | cycle 1 | cycle 2 | cycle 3 | cycle 4 | events that release energy |
|---|---|---|---|---|---|---|
| dx 3 | 30 | 63 % | 13 % | 7 % | 3 % | 25 of 26 |
| dx 1.5 | 30 | 43 % | 7 % | 7 % | 0 % | 17 of 17 |
