# Independent exact demag energy: Fourier sum over the aliases of the cube form factor (check_demag_alias.py)

Energy per box volume in K_d. Alias sum truncated at |G| ≤ Gm (4, 8, 16), extrapolated in 1/Gm. FD = DemagFieldPBC; Newell with 1, 2 (4) image shells and compensated outer surface charges. Run 2026-10-15 on the local CPU.

    cube f=1/8, m (111)          alias Gm (4, 8, 16): 0.0361075 0.0362724 0.0363625 -> extrap. 0.0364526 | FD 0.0381311 | Newell 0.0364583 0.0364583 | exact 0.0364583  [29 s]
    checkerboard (1,-1,0)        alias Gm (4, 8, 16): 0.4556193 0.4763894 0.4877866 -> extrap. 0.4991838 | FD 0.0000000 | Newell 0.5001063 0.5000383 | exact -  [34 s]
    B0b uniform (111)            alias Gm (4, 8, 16): 0.0754153 0.0757233 0.0758929 -> extrap. 0.0760626 | FD 0.0832963 | Newell 0.0760640 0.0760707 | exact 0.0760744  [243 s]
    B0b I3 s941 (box: FD/N4)     alias Gm (4, 8, 16): 0.0061062 0.0063603 0.0065003 -> extrap. 0.0066403 | FD 0.0008937 | Newell 0.0067204 0.0066748 0.0066578 | exact -  [335 s]
    B0b T3 s941 (box: FD/N4)     alias Gm (4, 8, 16): 0.0062331 0.0064874 0.0066275 -> extrap. 0.0067676 | FD 0.0008934 | Newell 0.0067926 0.0067822 0.0067788 | exact -  [328 s]
    [exited with code 0]
