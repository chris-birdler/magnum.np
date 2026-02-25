# %% [markdown]
# # Spin wave excited by a Rayleigh Surface Acoustic Wave (SAW)
# Reference: https://arxiv.org/abs/2507.23456

# %% [markdown]
# ## Run Simulation

# %%
import gc
import traceback
import torch
import magnumnp as mag
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from colorama import Fore, Style
def get_now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t"
MU_0 = mag.constants.mu_0
mag.Timer.enable()


# Define Simulation Parameters
save_dir = Path.cwd() / "data"

# External Field
# For the Rayleigh mode φ = n * 45° are the directions of maximum resonance
phi_vals = [45]  # [°] Angles in Degrees from propagation direction (x-axis)
# For the chosen material parameters (10nm Ni), the following field values correspond to (in order):
# main peak, secondary peak, off-resonance
H_vals = [-46, 46, 60]  # [mT] Field Strengths

mesh_n = [120, 2, 4]  # Number of Cells in each Direction
# Length of mesh in x direction gets determined by Frequency of SAW.
# It is set to be exactly 1 wavelength of the SAW.
mesh_size_y = 10e-9  # [m] Length of mesh in y direction: Irrelevant for ΔSij Calculation
mesh_size_z = 10e-9  # [m] Length of mesh in z direction

capping_thick = 5e-9  # [m] Thickness of capping layer: Does not get simulated directly,
                      # but influences at which depth the SAW interacts with the magnetic system

SAW_time_discret = 50  # [1] How many time discretizations there are during one SAW period:
                       # dt = 1 / (saw_t_d * f_saw)

simulation_time = 10e-9  # [s] Maximum simulation length
simulation_eps = 1e-5    # [1] Finishes the simulation early if the relative difference
                         # in energy absorption is less than this value [dB/dB]
simulation_abs = 1e-3    # [1] Finishes the simulation early if the absolute difference
                         # in energy absorption is less than this value [dB/mm]

log_every = 5  # [1] How often logs of magnetization and displacement are made.
               # If zero, no logs are made at all.


# SAW Parameters
SAW_frequency = 4.47e9  # [Hz] Frequency of the SAW
SAW_amplitude = 4e-12   # [m] Amplitude of SAW simulation. Has to be chosen sufficiently small
                        # to not enter non-linear regime and sufficiently large so as not to get
                        # numerical errors. Otherwise irrelevant for ΔSij Calculation.
reverse_direction = False  # Reverses direction of SAW


# Magnetic Material Parameters
# Nickel (Values are taken from https://doi.org/10.1103/PhysRevApplied.15.034046
# unless otherwise stated):
alpha = 69e-3          # [1] Gilbert Damping
A = 7.7e-12            # [J/m] Exchange Stiffness
Ms = 408e3             # [A/m] Saturation Magnetization
lambda_s = -14.2178e-6 # [1] Saturation Magnetostriction (Fit parameter)

H_surf = 93e3                      # [A/m] Surface Anisotropy (Fit parameter)
K_surf = H_surf * MU_0 * Ms / 2    # [J/m³] Surface Anisotropy (Fit parameter)

H_easy = 1.4e-3         # [T] In-Plane Easy Axis Anisotropy
H_easy_phi = 83.6       # [°] Angle of In-Plane Easy Axis Anisotropy
K_easy = H_easy * Ms / 2  # [J/m³] In-Plane Easy Axis Anisotropy

rho_mag = 8900   # [kg/m³] Density
E_mag = 218e9    # [Pa] Young's Modulus
nu_mag = 0.3     # [1] Poisson's Ratio
C_mag = mag.C_isotropic(E_mag, nu_mag)


# COMSOL Results for Phonon Energy
# [J] Potential + Kinetic + Electric Energy obtained from COMSOL simulation
# for A=9.791578e-21 m, f=3.109220097206838E9 Hz, 1µm length, 1m width
SAW_energy_comsol = (1.831145619991022E-29 + 1.841984428928603E-29 +
                     1.0838808936994005E-31)
SAW_amplitude_comsol = 9.791578e-21          # [m] Amplitude in COMSOL simulation
SAW_freq_comsol = 3.109220097206838E9        # [Hz] Frequency in COMSOL simulation
SAW_length_comsol = 1e-6                     # [m] Length of simulated SAW in COMSOL simulation
SAW_width_comsol = 1                         # [m] Width of simulated SAW in COMSOL simulation
SAW_ct = 3416.926985252044                   # [m/s] Obtained from COMSOL simulation (fit)
SAW_cl = 5545.0211171483315                  # [m/s] Obtained from COMSOL simulation (fit)

SAW_amp_factor = 0  # [m] If zero, AmpFactor will be calculated
SAW_energy = 0      # [J] If zero, Energy will be calculated from the parameters above




# SAW parameterization based on https://doi.org/10.1063/1.30437

def rayleigh_u(state, q, k_t, k_l, A, amp, omega, mesh_z, capping):
    """Calculate displacement field for Rayleigh wave."""
    nx, ny, nz = state.mesh.n
    u = torch.zeros((nx, ny, nz, 3))

    xx, xy, xz = state.mesh.SpatialCoordinate()
    # Conversion from magnum.np coordinate system to https://doi.org/10.1063/1.30437 coordinate system
    xz_down = mesh_z - xz + capping

    kt_exp = torch.exp(-k_t * xz_down)
    kl_exp = torch.exp(-k_l * xz_down)

    u[:, :, :, 0] = ((kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * kl_exp) *
                     k_t * torch.cos(xx * q - state.t * omega))

    u[:, :, :, 2] = ((kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * kl_exp) *
                     q * torch.sin(xx * q - state.t * omega))

    return u * (A / amp)


def rayleigh_v(state, q, k_t, k_l, A, amp, omega, mesh_z, capping):
    """Calculate velocity field for Rayleigh wave."""
    nx, ny, nz = state.mesh.n
    v = torch.zeros((nx, ny, nz, 3))

    xx, xy, xz = state.mesh.SpatialCoordinate()
    xz_down = mesh_z - xz + capping

    kt_exp = torch.exp(-k_t * xz_down)
    kl_exp = torch.exp(-k_l * xz_down)

    v[:, :, :, 0] = ((kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * kl_exp) *
                     k_t * omega * torch.sin(xx * q - state.t * omega))

    v[:, :, :, 2] = ((kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * kl_exp) *
                     q * (-omega) * torch.cos(xx * q - state.t * omega))

    return v * (A / amp)


def rayleigh_eps(state, q, k_t, k_l, A, amp, omega, mesh_z, capping):
    """Calculate strain tensor for Rayleigh wave."""
    nx, ny, nz = state.mesh.n
    eps = torch.zeros((nx, ny, nz, 6))

    xx, xy, xz = state.mesh.SpatialCoordinate()
    xz_down = mesh_z - xz + capping

    kt_exp = torch.exp(-k_t * xz_down)
    kl_exp = torch.exp(-k_l * xz_down)

    # eps_xx
    eps[:, :, :, 0] = ((kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * kl_exp) *
                       k_t * (-q) * torch.sin(xx * q - state.t * omega))
    # eps_yy
    eps[:, :, :, 1] = 0
    # eps_zz
    eps[:, :, :, 2] = (((k_t) * kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * (k_l) * kl_exp) *
                       q * torch.sin(xx * q - state.t * omega))
    # 2 * eps_yz
    eps[:, :, :, 3] = 0
    # 2 * eps_xz
    eps[:, :, :, 4] = (((k_t) * kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * (k_l) * kl_exp) *
                       k_t * torch.cos(xx * q - state.t * omega) +
                       (kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * kl_exp) *
                       q**2 * torch.cos(xx * q - state.t * omega))
    # 2 * eps_xy
    eps[:, :, :, 5] = 0

    return eps * (A / amp)


def rayleigh_deps_dt(state, q, k_t, k_l, A, amp, omega, mesh_z, capping):
    """Calculate time derivative of strain tensor for Rayleigh wave."""
    nx, ny, nz = state.mesh.n
    deps = torch.zeros((nx, ny, nz, 6))

    xx, xy, xz = state.mesh.SpatialCoordinate()
    xz_down = mesh_z - xz + capping

    kt_exp = torch.exp(-k_t * xz_down)
    kl_exp = torch.exp(-k_l * xz_down)

    # eps_xx
    deps[:, :, :, 0] = ((kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * kl_exp) *
                        k_t * (q * omega) * torch.cos(xx * q - state.t * omega))
    # eps_yy
    deps[:, :, :, 1] = 0
    # eps_zz
    deps[:, :, :, 2] = (((k_t) * kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * (k_l) * kl_exp) *
                        (-q * omega) * torch.cos(xx * q - state.t * omega))
    # 2 * eps_yz
    deps[:, :, :, 3] = 0
    # 2 * eps_xz
    deps[:, :, :, 4] = (((k_t) * kt_exp - ((2 * q**2) / (q**2 + k_t**2)) * (k_l) * kl_exp) *
                        k_t * omega * torch.sin(xx * q - state.t * omega) +
                        (kt_exp - ((2 * k_t * k_l) / (q**2 + k_t**2)) * kl_exp) *
                        q**2 * omega * torch.sin(xx * q - state.t * omega))
    # 2 * eps_xy
    deps[:, :, :, 5] = 0

    return deps * (A / amp)


def calc_xi(ct, cl):
    """
    Calculate xi parameter from elastic wave speeds.
    Based on "Theory of Elasticity" by Landau and Lifshitz; p121, eq. 24.13
    """
    b = ct / cl
    x = sp.symbols('x')

    term1 = -(8 - 48 * b**2)

    term2_eq = sp.Eq(45 * b**2 + 3 * sp.sqrt(3) * sp.sqrt(-64 * b**6 + 107 * b**4 - 62 * b**2 + 11) - 17, x**3)
    sol2 = sp.solve(term2_eq, x)
    term2 = [solu.evalf() for solu in sol2][0]  # Choose Real valued term

    xi_eq = sp.Eq(((term1 / (6 * term2)) + (2 / 3) * term2 + 8 / 3), x**2)
    solxi = sp.solve(xi_eq, x)
    xi = [solu.evalf() for solu in solxi][0]
    return float(abs(xi))  # Choose positive solution


# Calculate derived SAW parameters
SAW_dis_xi = calc_xi(SAW_ct, SAW_cl)
SAW_omega = SAW_frequency * 2.0 * np.pi
# Dispersion relation from "Theory of Elasticity" by Landau and Lifshitz; p121, eq. 24.12
SAW_k = SAW_omega / (SAW_ct * SAW_dis_xi)
SAW_wl = 2 * np.pi / SAW_k
SAW_c = SAW_wl * SAW_frequency  # = c_t * dis_xi
SAW_period = 1 / SAW_frequency

simulation_dt = 1 / (SAW_frequency * SAW_time_discret)

# "Theory of Elasticity" by Landau and Lifshitz; p120
SAW_kt = np.sqrt(SAW_k**2 - (SAW_omega / SAW_ct)**2)
SAW_kl = np.sqrt(SAW_k**2 - (SAW_omega / SAW_cl)**2)

if reverse_direction:
    SAW_k *= -1

# Assuming polycrystalline material:
lambda100 = lambda_s
lambda111 = lambda_s

mesh_dx = [SAW_wl / mesh_n[0], mesh_size_y / mesh_n[1], mesh_size_z / mesh_n[2]]
# Guarantees the simulation sees the region as a square, VERY IMPORTANT!
simulation_pbc = (1, int((3 * SAW_wl) / (2 * mesh_size_y)), 0)


# Files and Folders
log_name = str(save_dir / "Logs.dat")
results_name = str(save_dir / "Results.dat")

save_dir.mkdir(parents=True, exist_ok=True)

if not Path(log_name).exists():
    with open(log_name, "w") as lfile:
        lfile.write(f"{get_now()}Simulations started\n")
        lfile.write(f"phi_vals:\t{phi_vals}\nh_vals:\t{H_vals}\n\n")
else:
    with open(log_name, "a") as lfile:
        lfile.write(f"\n\n{get_now()}Simulations started\n")
        lfile.write(f"phi_vals:\t{phi_vals}\nh_vals:\t{H_vals}\n\n\n")

if not Path(results_name).exists():
    with open(results_name, "w") as ffile:
        ffile.write("phi [°]\t\th [mT]\t\tt [s]\t\t\t\tS_ij [dB/mm]\n\n")


# Calculate SAW Energy
# This first call to rayleigh_u is necessary, because in the Rayleigh formula from
# https://doi.org/10.1063/1.30437, A is only proportional to the amplitude,
# but is not the amplitude itself!
if SAW_amp_factor == 0:
    # mesh has to be sufficiently thick to capture maximum amplitude of SAW
    mesh_amp = mag.Mesh((mesh_n[0], mesh_n[1], 100 * mesh_n[2]), mesh_dx, pbc=(0, 0, 0))
    state_amp = mag.State(mesh_amp)
    SAW_amp_factor = float(torch.max(rayleigh_u(
        state_amp, q=SAW_k, k_t=SAW_kt, k_l=SAW_kl, A=1, amp=1,
        omega=SAW_omega, mesh_z=100 * mesh_size_z, capping=0
    )))
    del mesh_amp, state_amp
    gc.collect()
    torch.cuda.empty_cache()

if SAW_energy == 0:
    SAW_energy = (SAW_energy_comsol * (SAW_amplitude / SAW_amplitude_comsol)**2 *
                  (SAW_frequency / SAW_freq_comsol) * (SAW_wl / SAW_length_comsol) *
                  (mesh_size_y / SAW_width_comsol))

print(f"{get_now()}SAW Parameters:")
print(f"{get_now()}Wavelength:\t{SAW_wl * 1e9:.3f} nm")
print(f"{get_now()}Frequency:\t{SAW_frequency * 1e-9:.3f} GHz")
print(f"{get_now()}Amplitude:\t{SAW_amplitude:.3e} m")
print(f"{get_now()}AmpFactor:\t{SAW_amp_factor:.25e} m")
print(f"{get_now()}ETot of SAW:\t{SAW_energy:.25e} J")

with open(save_dir / "SAWEnergy.dat", 'w') as file:
    file.write(f"{get_now()}SAW wavelength:\t{SAW_wl * 1e9:.3f} nm\n")
    file.write(f"{get_now()}SAW frequency:\t{SAW_frequency * 1e-9:.3f} GHz\n")
    file.write(f"{get_now()}SAW Amplitude:\t{SAW_amplitude:.3e} m\n")
    file.write(f"{get_now()}SAW AmpFactor:\t{SAW_amp_factor:.25e} m\n")
    file.write(f"{get_now()}ETot of SAW:\t{SAW_energy:.25e} J\n")


# Do Actual Simulation

total_sims = len(phi_vals) * len(H_vals)
sim_counter = 1
flag_keyboard_interrupted = False


for phi_name in phi_vals:
    results = [] # Sim results for plotting in format: [(H, S_ij), ...]
    for h_name in H_vals:
        try:
            with open(log_name, "a") as lfile:
                lfile.write(f"{get_now()}Starting simulating phi = {phi_name:.3f} and h = {h_name:.3f}\n")
            
            logger = mag.FieldLogger(
                str(save_dir / f"Fields/fields_phi_{phi_name:.1f}_h_{h_name:.1f}.pvd"),
                ["m", "ud"],
                every=log_every
            )

            # Setting up mesh and material
            mesh = mag.Mesh(mesh_n, mesh_dx, pbc=simulation_pbc)
            state = mag.State(mesh)

            state.material = {
                "alpha": alpha,
                "A": A,
                "Ms": Ms,
                "C": C_mag,
                "lambda_100": lambda100,
                "lambda_111": lambda111,
                "Ku_surf": K_surf,
                "Ku_axis_surf": [0, 0, 1],
                "Ku_easy": K_easy,
                "Ku_axis_easy": [
                    np.cos(H_easy_phi * np.pi / 180),
                    np.sin(H_easy_phi * np.pi / 180),
                    0
                ],
                "rho": rho_mag,
            }

            # Defining Variables
            cell_size = state.mesh.dx[0] * state.mesh.dx[1] * state.mesh.dx[2]
            state.t = 0
            sim_steps = 0
            final_wave_count = 0
            sij_rel = 1e9
            sij_dif = 1e9
            sij_list = []
            t_list = []

            def print_info(state=state):
                print(f"{get_now()}Info for simulation {sim_counter} of {total_sims}:\tphi = {phi_name:.2f}°\th = {h_name:.2f}mT")
                print(f"{get_now()}Simulating in:\t\t\t{save_dir}")
                print(f"{get_now()}Simulating on:\t\t\t{state.device}")
                print(f"{get_now()}Current and previous Sij:\t{sij_list[-1]:.2f} dB\t{sij_list[-SAW_time_discret-1]:.2f} dB")
                print(f"{get_now()}Absolute difference and cutoff:\t{sij_dif:.1e} dB\t\t{simulation_abs:.1e} dB")
                print(f"{get_now()}Relative difference and cutoff:\t{1e2 * sij_rel:.5f}%\t\t{1e2 * simulation_eps:.5f}%")
                print(f"{get_now()}Current and max. time:\t\t{1e9 * state.t:.2f} ns\t\t{1e9 * simulation_time:.2f} ns")

            # Setting up magnetization
            phi_ext = phi_name * np.pi / 180
            h_ext = h_name * 1e-3 / MU_0

            state.m = state.Constant((
                np.sign(h_ext + 1e-9) * np.cos(phi_ext),
                np.sign(h_ext + 1e-9) * np.sin(phi_ext),
                0.
            ))

            # Setting up fields
            exchange = mag.ExchangeField()
            external = mag.ExternalField(h=[
                h_ext * np.cos(phi_ext),
                h_ext * np.sin(phi_ext),
                0
            ])
            sim_amp = SAW_amp_factor
            mag_el = mag.MagnetoElasticField(
                mechanical_strain=lambda state: rayleigh_eps(
                    state, q=SAW_k, k_t=SAW_kt, k_l=SAW_kl, A=SAW_amplitude,
                    amp=sim_amp, omega=SAW_omega, mesh_z=mesh_size_z,
                    capping=capping_thick
                )
            )
            uniax_easy = mag.UniaxialAnisotropyField(Ku="Ku_easy", Ku_axis="Ku_axis_easy")
            uniax_surface = mag.UniaxialAnisotropyField(Ku="Ku_surf", Ku_axis="Ku_axis_surf")
            demag = mag.DemagField(cache_dir=str(save_dir / "DemagKernel/"))

            # Relaxing Simulation
            print(f"{get_now()}Relaxing simulation {sim_counter} of {total_sims}:\tphi = {phi_name:.2f}°\th = {h_name:.2f}mT")
            print(f"{get_now()}Simulating in:\t{save_dir}")
            print(f"{get_now()}Simulating on:\t{state.device}")
            
            # Initially no SAW contribution for relaxation
            energy_contributions = [exchange, external, uniax_easy, uniax_surface, demag]
            llg = mag.LLGSolver(energy_contributions)
            llg.relax(state)
            state.t = 0

            # Starting Simulation
            print(f"{get_now()}Starting simulation {sim_counter} of {total_sims}:\tphi = {phi_name:.2f}°\th = {h_name:.2f}mT")
            print(f"{get_now()}Simulating in:\t{save_dir}")
            print(f"{get_now()}Simulating on:\t{state.device}")

            # Activating SAW contribution for simulation
            energy_contributions = [exchange, external, mag_el, uniax_easy, uniax_surface, demag]
            llg = mag.LLGSolver(energy_contributions)

            while state.t < simulation_time:
                llg.step(state, simulation_dt)
                t_list.append(float(state.t))

                if log_every != 0:
                    state.ud = rayleigh_u(
                        state, q=SAW_k, k_t=SAW_kt, k_l=SAW_kl, A=SAW_amplitude,
                        amp=sim_amp, omega=SAW_omega, mesh_z=mesh_size_z,
                        capping=capping_thick
                    )
                    logger << state
                sim_steps += 1

                # Calculate dE_Ph/dt = (d/dt e):C:e_m = (d/dt e):sig_m
                eps_m = mag.epsilon_m(state)
                sig_m = mag.sigma(state, eps_m)
                eps_v = rayleigh_deps_dt(
                    state, q=SAW_k, k_t=SAW_kt, k_l=SAW_kl, A=SAW_amplitude,
                    amp=sim_amp, omega=SAW_omega, mesh_z=mesh_size_z,
                    capping=capping_thick
                )
                dE_Ph = float(torch.sum(eps_v * sig_m * cell_size))

                # Calculate ΔSij:
                sij = (10 / np.log(10)) * (1e-3 / SAW_c) * (dE_Ph / SAW_energy)
                sij_list.append(sij)

                # Condition to end simulation early
                if sim_steps >= 2:
                    sij_dif = abs(sij_list[-1] - sij_list[-2])
                    sij_rel = abs(sij_dif / sij_list[-1])
                    if sij_rel < simulation_eps or sij_dif < simulation_abs:
                        final_wave_count += 1
                        # Only end when condition is valid over one period of the SAW
                        if final_wave_count > SAW_time_discret:
                            break
                    else:
                        final_wave_count = 0

                if sim_steps > SAW_time_discret and sim_steps % SAW_time_discret == 0:
                    print_info()

            sij_result = np.average(sij_list[-SAW_time_discret:])
            results.append((h_name, sij_result))

            # Saving Results
            if sim_steps > SAW_time_discret:
                print_info()
            print(f"{get_now()}Finished Simulation! Average Sij (Result):\t{sij_result:.3f} dB")

            with open(results_name, "a") as ffile:
                ffile.write(f"{phi_name:.2f}\t\t{h_name:.2f}\t\t{float(state.t):.25e}\t{sij_result:.25e}\n")

            with open(log_name, "a") as lfile:
                lfile.write(f"{get_now()}Finished simulating phi = {phi_name:.3f} and h = {h_name:.3f}\n")

            del mesh, state, logger
            gc.collect()
            torch.cuda.empty_cache()

            sim_counter += 1

        # Exception Handling
        except KeyboardInterrupt:
            flag_keyboard_interrupted = True
            with open(log_name, "a") as lfile:
                lfile.write(f"{get_now()}WARNING: Simulation aborted while working on phi = {phi_name:.3f} and h = {h_name:.3f}\n")
            break
        except Exception as err:
            with open(log_name, "a") as lfile:
                lfile.write(f"{get_now()}ERROR: Simulation Crashed while working on phi = {phi_name:.3f} and h = {h_name:.3f}\n")
            print(f"{Fore.RED}ERROR: Simulation Crashed while working on phi = {phi_name:.3f} and h = {h_name:.3f}")
            print(err)
            print(Style.RESET_ALL)
            
            error_dir = save_dir / "Errors"
            error_dir.mkdir(exist_ok=True)
            error_file = (error_dir / f"Error_phi_{phi_name:.1f}_h_{h_name:.1f}_{datetime.now().strftime('%Y_%m_%d___%H_%M_%S')}.txt")
            with open(error_file, 'w') as traceback_file:
                print(traceback.format_exc(), file=traceback_file)
            continue
    
    if flag_keyboard_interrupted:
        break


    # Plot results
    print(f"\n{get_now()}Generating plots...")

    # Extract data from results
    sim_h = np.array([r[0] for r in results])
    sim_sij = np.array([r[1] for r in results])

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Load expected results from file
    expected_results_file = Path(__file__).parent / "expected_results.dat"
    if expected_results_file.exists():
        with open(expected_results_file, 'r') as f:
            lines = f.readlines()[1:]  # Skip header
        
        # Parse all expected results for current phi (for the line)
        all_expected = []
        matching_expected = []
        for line in lines:
            if line.strip():
                parts = line.split()
                if len(parts) >= 4:
                    exp_phi = float(parts[0])
                    exp_h = float(parts[1])
                    exp_sij = float(parts[3])
                    
                    # Collect all expected results at this phi value
                    if abs(exp_phi - phi_name) < 0.01:
                        all_expected.append([exp_h, exp_sij])
                        
                        # Also check if this H value is in our simulation H_vals
                        if exp_h in H_vals:
                            matching_expected.append([exp_h, exp_sij])
        
        # Plot line connecting all expected results at this phi
        if all_expected:
            all_expected = np.array(all_expected)
            # Sort by H for proper line plotting
            sorted_indices = np.argsort(all_expected[:, 0])
            all_expected_sorted = all_expected[sorted_indices]
            ax.plot(all_expected_sorted[:, 0], all_expected_sorted[:, 1], 
                    color='#A23B72', linewidth=2.5, label='Expected', 
                    zorder=1, alpha=0.6)
        
        # Plot matching expected results as scatter points with low alpha
        if matching_expected:
            matching_expected = np.array(matching_expected)
            ax.scatter(matching_expected[:, 0], matching_expected[:, 1], 
                        s=400, c='#A23B72', marker='D', edgecolors='white', 
                        linewidths=2.5, zorder=2, alpha=0.4)
    
    # Plot simulation results as scatter points
    ax.scatter(sim_h, sim_sij, s=120, c='#2E86AB', marker='o', 
                edgecolors='white', linewidths=2, label='Simulation',
                zorder=3, alpha=0.9)

    # Styling
    ax.set_xlabel('External Field H [mT]', fontweight='bold', fontsize=12)
    ax.set_ylabel(r'$\Delta S_{ij}$ [dB/mm]', fontweight='bold', fontsize=12)
    ax.set_title(r'Magnetic Insertion Losses of Rayleigh SAW at $\phi = %.1f^\circ$' % phi_name, 
                fontweight='bold', fontsize=13)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', framealpha=0.9, fontsize=11)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

    plt.tight_layout()

    # Save figure
    plot_filename = save_dir / f"results_phi_{phi_name:.1f}.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"{get_now()}Plot saved to: {plot_filename}")

mag.Timer.print_report()