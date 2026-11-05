# PLAN 4.20 stage B: staircase geometry, 100 mT, seeds 941 ... 944

| arm | seed | done | relaxed conv. | iterations | Bloch relaxed | Bloch samples (cycles) | cycles | complete | w | P (W/cm^3) | mean dt (s) | run time (h) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C8 | 941 | True | True | 7000 | False | 0 [] | 9 | True | 8.0976e-06 | 217.5 | 6.52e-12 | 0.26 |
| C8 | 942 | True | True | 10000 | False | 0 [] | 7 | True | 1.1711e-05 | 314.5 | 6.7e-12 | 0.23 |
| C8 | 943 | True | True | 6000 | False | 0 [] | 7 | True | 4.9305e-06 | 132.4 | 6.66e-12 | 0.21 |
| C8 | 944 | True | True | 9000 | False | 0 [] | 7 | True | 2.1015e-05 | 564.4 | 6.67e-12 | 0.22 |
| D8 | 941 | True | True | 9500 | False | 42 [0, 3, 4, 5, 6, 7] | 8 | True | 7.9851e-05 | 2145 | 4.36e-12 | 0.54 |
| D8 | 942 | True | True | 11000 | False | 1 [3] | 8 | True | 3.108e-05 | 834.7 | 5.44e-12 | 0.47 |
| D8 | 943 | True | True | 6000 | False | 0 [] | 8 | True | 6.5757e-06 | 176.6 | 6.7e-12 | 0.36 |
| D8 | 944 | True | True | 9500 | False | 24 [2, 3, 4, 5, 6, 7] | 8 | True | 4.7848e-05 | 1285 | 5.11e-12 | 0.47 |
| D8V1 | 941 | no run |  |  |  |  |  |  |  |  |  |  |
| D8V1 | 942 | no run |  |  |  |  |  |  |  |  |  |  |
| D8V1 | 943 | no run |  |  |  |  |  |  |  |  |  |  |
| D8V1 | 944 | no run |  |  |  |  |  |  |  |  |  |  |
| D8V2 | 941 | True | True | 10500 | False | 0 [] | 8 | True | 7.0774e-05 | 1901 | 5.93e-12 | 0.50 |
| D8V2 | 942 | True | True | 3500 | False | 0 [] | 8 | True | 3.1875e-05 | 856.1 | 6.42e-12 | 0.41 |
| D8V2 | 943 | True | True | 13000 | False | 0 [] | 8 | True | 6.9203e-06 | 185.9 | 6.68e-12 | 0.47 |
| D8V2 | 944 | True | True | 3000 | False | 0 [] | 8 | True | 4.6694e-05 | 1254 | 6.33e-12 | 0.42 |
| F8 | 941 | True | True | 11000 | False | 0 [] | 16 | True | - | - | 2.11e-12 | 2.17 |
| F8 | 942 | True | True | 28000 | False | 0 [] | 8 | True | 5.2403e-05 | 1407 | 2.11e-12 | 1.68 |
| F8 | 943 | True | True | 18500 | False | 0 [] | 10 | True | 3.753e-05 | 1008 | 2.12e-12 | 1.60 |
| F8 | 944 | True | True | 19000 | False | 0 [] | 8 | True | 4.3066e-05 | 1157 | 2.11e-12 | 1.42 |

- B-0: F8 has values in 3 seeds and Bloch points in 0 seeds: valid reference.
- B-a: D8 has Bloch samples in seeds [941, 942, 944]: B-a valid.
  - V2: ran its cycles in [941, 942, 943, 944], Bloch in none: passes

## B-b: overlap with the f-weighted F8 reference

| arm | seed | cut-cell overlap per particle | mean | cut cells > 45 deg | full cells next to the cut layer > 45 deg | cut cells left out |
|---|---|---|---|---|---|---|
| C8 | 941 | 0.8488 0.5524 0.6680 0.9871 | 0.7641 | 4129 | 2249 | 0 |
| C8 | 942 | 0.9115 0.9841 0.4430 0.9811 | 0.8299 | 2917 | 1628 | 0 |
| C8 | 943 | 0.9577 0.8542 0.8644 0.2970 | 0.7433 | 3750 | 2114 | 0 |
| C8 | 944 | 0.9255 0.8183 0.8583 0.9551 | 0.8893 | 1148 | 740 | 0 |
| D8 | 941 | 0.9061 0.6048 0.7373 0.9977 | 0.8115 | 3098 | 1668 | 0 |
| D8 | 942 | 0.9494 0.9968 0.5078 0.9968 | 0.8627 | 2565 | 1390 | 0 |
| D8 | 943 | 0.9817 0.9095 0.9226 0.4006 | 0.8036 | 3083 | 1720 | 0 |
| D8 | 944 | 0.9583 0.8665 0.9048 0.9769 | 0.9267 | 702 | 446 | 0 |
| D8V2 | 941 | 0.8817 0.6062 0.7200 0.9971 | 0.8013 | 3266 | 1761 | 0 |
| D8V2 | 942 | 0.9505 0.9962 0.5029 0.9932 | 0.8607 | 2576 | 1398 | 0 |
| D8V2 | 943 | 0.9789 0.9123 0.9251 0.3499 | 0.7916 | 3210 | 1778 | 0 |
| D8V2 | 944 | 0.9587 0.8619 0.8809 0.9735 | 0.9187 | 807 | 484 | 0 |

## B-c: factors on P (exp of the mean of ln over the seeds, 95 % t-interval)

- C8/F8: 0.243 [0.047, 1.252], n 3, SE of ln 0.381
- C8/D8: 0.335 [0.087, 1.294], n 4, SE of ln 0.425
- D8/F8: 0.487 [0.047, 5.020], n 3, SE of ln 0.542
- D8V2/F8: 0.495 [0.053, 4.673], n 3, SE of ln 0.522
- D8V2/D8: 0.983 [0.871, 1.109], n 4, SE of ln 0.038
- F8/D8: 2.054 [0.199, 21.172], n 3, SE of ln 0.542

**Selection for stage C: V2** (only V2 passes)
