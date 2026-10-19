# Stage A of PLAN.md: static test of the near-field correction (V1, Donahue & McMichael 2007)

## Geometry d300_phi0.62, N 152, q 2 (static_A_d300_phi0.62_q2.json)

R0: passed. Voxel cells 2179360, magnetic V1 cells 2225056, modified 98304, new 45696; φ_vox 0.62058, φ_f 0.61976.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.7e-10, voxel -3.1e-10; uniform_110 ref -2.2e-10, voxel +2.8e-10; uniform_111 ref -6.8e-11, voxel +8.8e-10.

- uniform_100: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +5.5e-11; P2 +5.5e-11; P3 +5.5e-11
- uniform_110: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +5.0e-10; P2 +5.0e-10; P3 +5.0e-10
- uniform_111: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +9.4e-10; P2 +9.4e-10; P3 +9.4e-10

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K2_0_nx_pd300 | -1.042e-02 | 4.719e-04 | 0.045 | -2.113e-04 | 3.325e-02 | 1.364e-02 |
| K2_1_nx_pd300 | -1.040e-02 | 4.716e-04 | 0.045 | -2.120e-04 | 3.330e-02 | 1.350e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 1 | 0.948 / 0.0429 / 0.169; 0.948 / 0.0430 / 0.168 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.594 / 0.0269 / 0.096; 0.594 / 0.0269 / 0.096 | E_K, e_E, e_h | no | not computed |
| P | 3 | 0.432 / 0.0196 / 0.067; 0.431 / 0.0196 / 0.067 | E_K, e_E, e_h | no | not computed |
| L | 1 | 0.952 / 0.0432 / 0.194; 0.953 / 0.0432 / 0.194 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.602 / 0.0273 / 0.133; 0.602 / 0.0273 / 0.133 | E_K, e_E, e_h | no | not computed |
| L | 3 | 0.442 / 0.0200 / 0.113; 0.442 / 0.0200 / 0.112 | E_K, e_E, e_h | no | not computed |

## Geometry d300_phi0.55, N 156, q 2 (static_A_d300_phi0.55_q2.json)

R0: passed. Voxel cells 2089600, magnetic V1 cells 2135776, modified 93792, new 46176; φ_vox 0.55041, φ_f 0.55010.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.0e-10, voxel -2.4e-10; uniform_110 ref -1.5e-10, voxel +3.8e-10; uniform_111 ref +9.8e-12, voxel +1.0e-09.

- uniform_100: E_ref 0.082497, E_vox 0.082486, E_frac 0.080952; window error of P (relative): P1 +6.1e-11; P2 +6.1e-11; P3 +6.1e-11
- uniform_110: E_ref 0.082497, E_vox 0.082486, E_frac 0.080952; window error of P (relative): P1 +5.3e-10; P2 +5.3e-10; P3 +5.3e-10
- uniform_111: E_ref 0.082497, E_vox 0.082486, E_frac 0.080952; window error of P (relative): P1 +9.9e-10; P2 +9.9e-10; P3 +9.9e-10

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K4_1_nx_pd300 | -1.847e-03 | 4.209e-04 | 0.228 | -1.894e-04 | 3.382e-02 | 1.264e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 1 | 0.939 / 0.2140 / 0.161 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.588 / 0.1340 / 0.091 | E_K, e_E, e_h | no | not computed |
| P | 3 | 0.428 / 0.0974 / 0.064 | E_K, e_E, e_h | no | not computed |
| L | 1 | 0.945 / 0.2153 / 0.184 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.597 / 0.1361 / 0.125 | E_K, e_E, e_h | no | not computed |
| L | 3 | 0.440 / 0.1001 / 0.105 | E_K, e_E, e_h | no | not computed |

