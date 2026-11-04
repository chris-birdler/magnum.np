**PLAN 4.18 step 3** (seeds 921, 923, 925; arms C, D, Dk at 100 mT)

| arm | seed | relax iterations | samples with a Bloch point, cycles >= 1 | cycles >= 2 | mean dt (ps), cycles >= 1 | w |
|---|---|---|---|---|---|---|
| C | 921 | 11000 | 0 of 1542 | 0 | 5.87 | 1.382e-05 |
| C | 923 | 15000 | 0 of 1542 | 0 | 5.85 | 1.525e-05 |
| C | 925 | 14000 | 0 of 1542 | 0 | 5.85 | 1.242e-05 |
| D | 921 | 25500 | 48 of 1799 | 37 | 2.91 | 2.213e-05 |
| D | 923 | 15500 | 5 of 1542 | 4 | 3.74 | 1.777e-05 |
| D | 925 | 31000 | 13 of 1542 | 12 | 3.41 | 1.827e-05 |
| Dk | 921 | 12500 | 0 of 1542 | 0 | 2.39 | 2.547e-05 |
| Dk | 923 | 6000 | 0 of 1542 | 0 | 2.39 | 2.298e-05 |
| Dk | 925 | 26500 | 0 of 1542 | 0 | 2.41 | 1.654e-05 |

3a (D): 66 samples with a Bloch point in the cycles >= 1; share with a cut cell in the cube 1.000; share with a large-angle pair in the cube 1.000; holds: True
3b: seeds with a Bloch point of D in the cycles >= 2: {921: True, 923: True, 925: True}; fails 0; holds: True
3c: seeds that take part: {921: False, 923: False, 925: False}; holds: False
3f D/C: factor 1.400 [0.930, 2.106] (exp of the mean ln, SE of ln 0.095, 95 % t-interval, n 3); per seed {921: 1.601, 923: 1.165, 925: 1.47}
3f Dk/C: factor 1.546 [1.029, 2.323] (exp of the mean ln, SE of ln 0.095, 95 % t-interval, n 3); per seed {921: 1.842, 923: 1.506, 925: 1.332}
3f Dk/D: factor 1.105 [0.704, 1.733] (exp of the mean ln, SE of ln 0.105, 95 % t-interval, n 3); per seed {921: 1.151, 923: 1.293, 925: 0.906}
3g C: pooled {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}; per seed {921: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}, 923: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}, 925: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}}
3g D: pooled {'thin': 66, 'full': 0, 'cut': 0, 'rows': 66, 'events': 47}; per seed {921: {'thin': 48, 'full': 0, 'cut': 0, 'rows': 48, 'events': 31}, 923: {'thin': 5, 'full': 0, 'cut': 0, 'rows': 5, 'events': 5}, 925: {'thin': 13, 'full': 0, 'cut': 0, 'rows': 13, 'events': 11}}
3g Dk: pooled {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}; per seed {921: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}, 923: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}, 925: {'thin': 0, 'full': 0, 'cut': 0, 'rows': 0, 'events': 0}}
3g place of D: place 1 (seeds with at least 5 rows: [921, 923, 925])

Verdict of step 3: not determined
