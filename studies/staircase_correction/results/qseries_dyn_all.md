# Dynamic series in q: readings of PLAN 9 (all runs, also the runs with an event; information only)

| arm | runs done | runs with an event | runs without an event record | geometric mean of P (W/cm³), runs used |
|---|---|---|---|---|
| T3 | 12 | 3 | 0 | 1.01 (n 12) |
| T3Pc | 12 | 4 | 0 | 1.34 (n 12) |
| T3Pc4 | 12 | 6 | 0 | 1.15 (n 12) |
| X3q4 | 5 | 2 | 7 | 0.747 (n 5) |
| T3Pc8 | 0 | 0 | 12 | – |

## Overlap of the relaxed states (init.pt), per particle, common magnetic cells

| pair | seed | particle 0 / 1 / 2 / 3 |
|---|---|---|
| T3Pc4 − T3Pc | 941 | 0.9950 / 0.9982 / 0.9978 / 0.9990 |
| T3Pc4 − T3Pc | 942 | 0.9962 / 0.9992 / 0.9991 / 0.9985 |
| T3Pc4 − T3Pc | 943 | 0.9985 / 0.9992 / 0.9975 / 0.9959 |
| T3Pc4 − T3Pc | 944 | 0.9976 / 0.9957 / 0.9994 / 0.9987 |
| T3Pc4 − T3Pc | 945 | 0.9992 / 0.9994 / 0.9983 / 0.9979 |
| T3Pc4 − T3Pc | 946 | 0.9960 / 0.9992 / 0.9983 / 0.9964 |
| T3Pc4 − T3Pc | 947 | 0.9981 / 0.9996 / 0.9972 / 0.9990 |
| T3Pc4 − T3Pc | 948 | 0.9967 / 0.9987 / 0.9996 / 0.9970 |
| T3Pc4 − T3Pc | 949 | 0.9988 / 0.9993 / 0.9951 / 0.9996 |
| T3Pc4 − T3Pc | 950 | 0.9962 / 0.9907 / 0.9987 / 0.9990 |
| T3Pc4 − T3Pc | 951 | 0.9969 / 0.9976 / 0.9996 / 0.9995 |
| T3Pc4 − T3Pc | 952 | 0.9990 / 0.9983 / 0.9948 / 0.9951 |
| T3Pc4 − X3q4 | 941 | 1.0000 / 1.0000 / 1.0000 / 1.0000 |
| T3Pc4 − X3q4 | 942 | 1.0000 / 1.0000 / 1.0000 / 1.0000 |
| T3Pc4 − X3q4 | 943 | 1.0000 / 1.0000 / 1.0000 / 1.0000 |
| T3Pc4 − X3q4 | 944 | 1.0000 / 1.0000 / 1.0000 / 1.0000 |
| T3Pc4 − X3q4 | 945 | 1.0000 / 1.0000 / 1.0000 / 1.0000 |

## Effects (ln P; mean ± SE; 95 % t-interval; band ±0.15)

| effect | meaning | per seed 941 / 942 / 943 / 944 / 945 / 946 / 947 / 948 / 949 / 950 / 951 / 952 | mean ± SE (n) | 95 % interval | reading |
|---|---|---|---|---|---|
| N3(q 4) | correction against exact field at q 4: ln P(T3Pc4) − ln P(X3q4) | -0.0006 / +0.0238 / +0.0015 / +0.0028 / -0.0005 / – / – / – / – / – / – / – | +0.0054 ± 0.0046 (SE, n 5) | -0.0075 … +0.0183 | within ±0.15 |
| Δ_24 | step q 2 → 4: ln P(T3Pc4) − ln P(T3Pc) | -0.4540 / -0.2959 / -0.0521 / -0.6537 / -0.0379 / -0.0617 / -0.1122 / -0.0383 / +0.8553 / -1.3294 / -0.0697 / +0.4291 | -0.1517 ± 0.1546 (SE, n 12) | -0.4920 … +0.1885 | not determined (information) |
| 2 Δ_24 | step q 2 → smooth spheres (first order in 1/q), phase A | -0.9081 / -0.5919 / -0.1042 / -1.3075 / -0.0757 / -0.1234 / -0.2243 / -0.0765 / +1.7105 / -2.6588 / -0.1394 / +0.8581 | -0.3034 ± 0.3092 (SE, n 12) | -0.9839 … +0.3771 | not determined |
| Δ_c(4) | correction at q 4 against Coarse: ln P(T3Pc4) − ln P(T3) | +0.0743 / +0.0074 / -0.0429 / +0.0418 / -0.0504 / -2.3581 / +0.0816 / -0.0507 / +2.4850 / +0.0515 / +0.0017 / +1.2394 | +0.1234 ± 0.3152 (SE, n 12) | -0.5703 … +0.8171 | not determined (information) |
| Δ_48 | step q 4 → 8: ln P(T3Pc8) − ln P(T3Pc4) | – / – / – / – / – / – / – / – / – / – / – / – | n 0 | – | not determined (n 0) (information) |
| Δ_24 + 2 Δ_48 | step q 2 → smooth spheres (first order in 1/q), phase B | – / – / – / – / – / – / – / – / – / – / – / – | n 0 | – | not determined (n 0) |
| Δ_48 − Δ_24 / 2 | check of first order in 1/q | – / – / – / – / – / – / – / – / – / – / – / – | n 0 | – | not determined (n 0) |

Start rule of phase B (PLAN 9.3): not applied here (it uses the event-free pairs).
