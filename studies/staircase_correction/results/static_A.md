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

## Geometry d300_phi0.62, N 152, q 2 (static_A_d300_phi0.62_q2_p456.json)

R0: passed. Voxel cells 2179360, magnetic V1 cells 2225056, modified 98304, new 45696; φ_vox 0.62058, φ_f 0.61976.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.7e-10, voxel -3.1e-10; uniform_110 ref -2.2e-10, voxel +2.8e-10; uniform_111 ref -6.8e-11, voxel +8.8e-10.

- uniform_100: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P4 +5.5e-11; P5 +5.5e-11; P6 +5.5e-11
- uniform_110: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P4 +5.0e-10; P5 +5.0e-10; P6 +5.0e-10
- uniform_111: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P4 +9.4e-10; P5 +9.4e-10; P6 +9.4e-10

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K2_0_nx_pd300 | -1.042e-02 | 4.719e-04 | 0.045 | -2.113e-04 | 3.325e-02 | 1.364e-02 |
| K2_1_nx_pd300 | -1.040e-02 | 4.716e-04 | 0.045 | -2.120e-04 | 3.330e-02 | 1.350e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 4 | 0.339 / 0.0154 / 0.052; 0.338 / 0.0154 / 0.052 | E_K, e_E, e_h | no | not computed |
| P | 5 | 0.278 / 0.0126 / 0.043; 0.278 / 0.0126 / 0.042 | E_K, e_E, e_h | no | not computed |
| P | 6 | 0.236 / 0.0107 / 0.036; 0.236 / 0.0107 / 0.036 | E_K, e_E, e_h | no | not computed |
| L | 4 | 0.352 / 0.0159 / 0.103; 0.351 / 0.0159 / 0.102 | E_K, e_E, e_h | no | not computed |
| L | 5 | 0.293 / 0.0133 / 0.097; 0.293 / 0.0133 / 0.096 | E_K, e_E, e_h | no | not computed |
| L | 6 | 0.252 / 0.0114 / 0.094; 0.251 / 0.0114 / 0.093 | E_K, e_E, e_h | no | not computed |

## Geometry d300_phi0.62, N 152, q 2 (static_A_d300_phi0.62_q2_c.json)

R0: passed. Voxel cells 2179360, magnetic V1 cells 2225056, modified 98304, new 45696; φ_vox 0.62058, φ_f 0.61976.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.7e-10, voxel -3.1e-10; uniform_110 ref -2.2e-10, voxel +2.8e-10; uniform_111 ref -6.8e-11, voxel +8.8e-10.

- uniform_100: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +5.5e-11; P2 +5.5e-11; P3 +5.5e-11
- uniform_110: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +5.0e-10; P2 +5.0e-10; P3 +5.0e-10
- uniform_111: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P1 +9.5e-10; P2 +9.5e-10; P3 +9.5e-10

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K2_0_nx_pd300 | -1.042e-02 | 4.719e-04 | 0.045 | -2.113e-04 | 3.325e-02 | 1.364e-02 |
| K2_1_nx_pd300 | -1.040e-02 | 4.716e-04 | 0.045 | -2.120e-04 | 3.330e-02 | 1.350e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 1 | 0.193 / 0.0088 / 0.047; 0.193 / 0.0087 / 0.047 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.049 / 0.0022 / 0.011; 0.049 / 0.0022 / 0.011 | E_K, e_E, e_h | yes | not computed |
| P | 3 | 0.018 / 0.0008 / 0.004; 0.018 / 0.0008 / 0.004 | E_K, e_E, e_h | yes | not computed |
| L | 1 | 0.189 / 0.0086 / 0.080; 0.188 / 0.0085 / 0.080 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.045 / 0.0020 / 0.062; 0.044 / 0.0020 / 0.062 | E_K, e_E, e_h | yes | not computed |
| L | 3 | 0.014 / 0.0006 / 0.061; 0.013 / 0.0006 / 0.061 | E_K, e_E, e_h | yes | not computed |

## Geometry d300_phi0.55, N 156, q 2 (static_A_d300_phi0.55_q2_c.json)

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
| P | 1 | 0.196 / 0.0447 / 0.046 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.050 / 0.0115 / 0.010 | E_K, e_E, e_h | yes | not computed |
| P | 3 | 0.019 / 0.0043 / 0.004 | E_K, e_E, e_h | yes | not computed |
| L | 1 | 0.194 / 0.0442 / 0.076 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.048 / 0.0110 / 0.058 | E_K, e_E, e_h | yes | not computed |
| L | 3 | 0.016 / 0.0037 / 0.056 | E_K, e_E, e_h | yes | not computed |

## Geometry d300_phi0.62, N 152, q 2 (static_A_d300_phi0.62_q2_k1_c.json)

R0: passed. Voxel cells 2179360, magnetic V1 cells 2225056, modified 98304, new 45696; φ_vox 0.62058, φ_f 0.61976.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.7e-10, voxel -3.1e-10; uniform_110 ref -1.5e-06, voxel -6.1e-06; uniform_111 ref -3.1e-06, voxel -1.2e-05.

- uniform_100: E_ref 0.078553, E_vox 0.078487, E_frac 0.076804; window error of P (relative): P2 +5.5e-11; P3 +5.5e-11
- uniform_110: E_ref 0.078553, E_vox 0.078486, E_frac 0.076803; window error of P (relative): P2 -4.6e-06; P3 -4.6e-06
- uniform_111: E_ref 0.078553, E_vox 0.078486, E_frac 0.076803; window error of P (relative): P2 -9.2e-06; P3 -9.2e-06

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K2_0_nx_pd300 | -1.042e-02 | 4.719e-04 | 0.045 | -2.113e-04 | 3.325e-02 | 1.364e-02 |
| K2_1_nx_pd300 | -1.040e-02 | 4.716e-04 | 0.045 | -2.120e-04 | 3.330e-02 | 1.350e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 2 | 0.049 / 0.0022 / 0.011; 0.049 / 0.0022 / 0.011 | E_K, e_E, e_h | yes | not computed |
| P | 3 | 0.018 / 0.0008 / 0.004; 0.018 / 0.0008 / 0.004 | E_K, e_E, e_h | yes | not computed |
| L | 2 | 0.045 / 0.0020 / 0.062; 0.044 / 0.0020 / 0.062 | E_K, e_E, e_h | yes | not computed |
| L | 3 | 0.014 / 0.0006 / 0.061; 0.013 / 0.0006 / 0.061 | E_K, e_E, e_h | yes | not computed |

## Geometry d300_phi0.62, N 152, q 4 (static_A_d300_phi0.62_q4_k1_c.json)

R0: passed. Voxel cells 2179360, magnetic V1 cells 2251392, modified 144896, new 72032; φ_vox 0.62058, φ_f 0.62001.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -4.2e-10, voxel -3.1e-10; uniform_110 ref -3.8e-07, voxel -6.1e-06; uniform_111 ref -7.7e-07, voxel -1.2e-05.

- uniform_100: E_ref 0.078532, E_vox 0.078487, E_frac 0.076571; window error of P (relative): P1 +1.1e-10; P2 +1.1e-10; P3 +1.1e-10
- uniform_110: E_ref 0.078532, E_vox 0.078486, E_frac 0.076570; window error of P (relative): P1 -5.7e-06; P2 -5.7e-06; P3 -5.7e-06
- uniform_111: E_ref 0.078532, E_vox 0.078486, E_frac 0.076570; window error of P (relative): P1 -1.1e-05; P2 -1.1e-05; P3 -1.1e-05

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K2_0_nx_pd300 | -1.042e-02 | 7.957e-04 | 0.076 | 3.783e-05 | 3.529e-02 | 1.082e-02 |
| K2_1_nx_pd300 | -1.040e-02 | 7.962e-04 | 0.077 | 3.763e-05 | 3.532e-02 | 1.065e-02 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 1 | 0.152 / 0.0116 / 0.045; 0.152 / 0.0116 / 0.045 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.039 / 0.0030 / 0.011; 0.039 / 0.0030 / 0.011 | E_K, e_E, e_h | yes | not computed |
| P | 3 | 0.014 / 0.0011 / 0.004; 0.014 / 0.0011 / 0.004 | E_K, e_E, e_h | yes | not computed |
| L | 1 | 0.149 / 0.0114 / 0.078; 0.149 / 0.0114 / 0.078 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.036 / 0.0027 / 0.059; 0.035 / 0.0027 / 0.059 | E_K, e_E, e_h | yes | not computed |
| L | 3 | 0.011 / 0.0009 / 0.057; 0.011 / 0.0008 / 0.057 | E_K, e_E, e_h | yes | not computed |

## Geometry d300_phi0.55, N 156, q 4 (static_A_d300_phi0.55_q4_k1_c.json)

R0: passed. Voxel cells 2089600, magnetic V1 cells 2159200, modified 142208, new 69600; φ_vox 0.55041, φ_f 0.55003.

R1: passed. Relative deviation from φ(1 − φ)/3: uniform_100 ref -3.5e-10, voxel -2.4e-10; uniform_110 ref -2.8e-07, voxel -4.4e-06; uniform_111 ref -5.5e-07, voxel -8.9e-06.

- uniform_100: E_ref 0.082499, E_vox 0.082486, E_frac 0.080728; window error of P (relative): P1 +1.1e-10; P2 +1.1e-10; P3 +1.1e-10
- uniform_110: E_ref 0.082499, E_vox 0.082486, E_frac 0.080728; window error of P (relative): P1 -4.2e-06; P2 -4.2e-06; P3 -4.2e-06
- uniform_111: E_ref 0.082499, E_vox 0.082485, E_frac 0.080728; window error of P (relative): P1 -8.3e-06; P2 -8.3e-06; P3 -8.3e-06

| state | E_K | S = E_vox − E_ref | S / |E_K| | E_frac − E_ref | RMS h_vox − h_ref (layer) | RMS h_frac − h_ref |
|---|---|---|---|---|---|---|
| K4_1_nx_pd300 | -1.847e-03 | 7.104e-04 | 0.385 | 3.085e-05 | 3.521e-02 | 9.728e-03 |

| form | p | per state: e_E / |ΔE|/|E_K| / e_h | criteria in use | R2 in all states | R3 (λ_min / λ_max) |
|---|---|---|---|---|---|
| P | 1 | 0.157 / 0.0605 / 0.045 | E_K, e_E, e_h | no | not computed |
| P | 2 | 0.040 / 0.0155 / 0.011 | E_K, e_E, e_h | yes | not computed |
| P | 3 | 0.015 / 0.0058 / 0.004 | E_K, e_E, e_h | yes | not computed |
| L | 1 | 0.156 / 0.0599 / 0.077 | E_K, e_E, e_h | no | not computed |
| L | 2 | 0.039 / 0.0149 / 0.057 | E_K, e_E, e_h | yes | not computed |
| L | 3 | 0.013 / 0.0051 / 0.055 | E_K, e_E, e_h | yes | not computed |

