#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import glob
import os
import subprocess
import sys
from scipy.optimize import curve_fit

# -----------------------------
# Read MONKES database
# -----------------------------
# This function reads the output of MONKES after it has finished running.
# It extracts geometric quantities and monoenergetic transport coefficients,
# fits models for D11, and appends all information to monkes.dk.
def read_monkes_database(params):
    import glob

    print("\n--- Interpolating 7 surfaces for MONKES ---")
    s = params['s']
    iota_arr = params['iota_arr']
    R00 = params['R00']
    minorradiusVMEC = params['minorradiusVMEC']
    Nperiods = params['Nperiods']

    with open("monkes.dk", "a", encoding="utf-8") as f:
        # Loop through the 7 MONKES surfaces
        for surf_idx in range(1, 8):
            # Compute radial coordinate and interpolated quantities
            s_val = ((surf_idx - 0.5) / 7.0)**2
            if surf_idx == 1:
                s_val *= 3.0
            r_val    = np.sqrt(s_val) * minorradiusVMEC
            iota_val = np.interp(s_val, s, iota_arr)
            B_val_1T = 1.0
            R_val    = np.interp(s_val, s, R00)

            print(f"\nSurface {surf_idx}: s_val={s_val:.6e}, r_val={r_val:.6e}, "
                  f"iota_val={iota_val:.6e}, B_val_1T={B_val_1T:.6e}, R_val={R_val:.6e}")

            folder_s = f"s_{s_val:.3f}"

            # Initialize values for geometric computations
            b2av = 1.0
            xkn = 0.0
            ft = 0.0
            b10 = 0.0
            bmn_dict = {}

            # Path to MONKES-generated magnetic configuration
            config_path = os.path.join(folder_s, "cmul_1.0e+00", "monkes_Magnetic_configuration.dat")
            if os.path.exists(config_path):
                print(f"Reading MONKES config: {config_path}")
                with open(config_path, "r") as cfg:
                    for line in cfg:
                        cols = line.split()
                        if len(cols) != 4:
                            continue
                        try:
                            # Read Fourier coefficients
                            m = int(cols[0])
                            n = int(cols[1])
                            ratio = float(cols[3])
                            bmn_dict[(m,n)] = ratio
                            if m == 1 and n == 0:
                                b10 = ratio
                            if m == 0 and n == 0:
                                B_val = float(cols[2])
                        except ValueError:
                            continue

                # Construct a grid to reconstruct B(θ,ζ)
                theta = np.linspace(0, 2*np.pi, 256)
                zeta = np.linspace(0, 2*np.pi, 256)
                THETA, ZETA = np.meshgrid(theta, zeta, indexing='ij')

                # Sum Fourier modes to build field strength on the surface
                b = np.zeros_like(THETA)
                for (m,n), bmn in bmn_dict.items():
                    b += bmn * np.cos(m*THETA - Nperiods*n*ZETA)

                # Compute <b^2> using the harmonic average
                total_points = b.size
                sum1_over_b2 = np.sum(1.0 / b**2)
                b2av = total_points / sum1_over_b2

                # Compute the trapped fraction ft
                B_max = B_val * np.max(b)
                Nl = 4000
                lambda_max = 1.0 / B_max
                lambdas = np.linspace(0, lambda_max, Nl)
                inv_B2 = 1.0 / (B_val * b)**2
                den_sum = np.sum(inv_B2)
                sqrt1mlamB_av = np.zeros_like(lambdas)

                # Integrate over pitch-angle parameter λ
                for i, lam in enumerate(lambdas):
                    arg = 1 - lam * B_val * b
                    arg = np.maximum(arg, 0.0)
                    num = np.sum(np.sqrt(arg) * inv_B2)
                    sqrt1mlamB_av[i] = num / den_sum

                integrand = lambdas / sqrt1mlamB_av
                ft = 1 - (3 * B_val**2 * b2av / 4) * np.trapz(integrand, lambdas)

                # Compute normalized mirror term xkn
                xkn = np.abs(b10) / (r_val / R_val)

                # Write basic geometry invariants for this surface
                f.write(f"{r_val:14.6e} {R_val:14.6e} {B_val_1T:14.6e} {iota_val:14.6e} "
                        f"{xkn:14.6f} {ft:14.6f} {b2av:14.6f} r,R,B,io,xkn,ft,<b^2>\n")

                # Read all MONKES monoenergetic databases for this surface
                cmul_all, efield_all, Ntz_all, Nxi_all, D11_all, D31_all, D33_all = [],[],[],[],[],[],[]
                if os.path.exists(folder_s) and os.path.isdir(folder_s):
                    for sub_folder in sorted(glob.glob(os.path.join(folder_s,"*"))):
                        db_files = glob.glob(os.path.join(sub_folder,"monkes_Monoenergetic_Database.dat"))
                        for db_file in db_files:
                            print(f"Reading MONKES database: {db_file}")
                            with open(db_file,"r") as db:
                                lines_db = db.readlines()[1:] # skip header
                                for line_db in lines_db:
                                    cols = line_db.split()
                                    if len(cols)<9:
                                        continue
                                    # Collect transport coefficients
                                    cmul_all.append(float(cols[0]))
                                    efield_all.append(float(cols[1]))
                                    Ntz_all.append(float(cols[2])*float(cols[3]))
                                    Nxi_all.append(float(cols[4]))
                                    D11_all.append(float(cols[5]))
                                    D31_all.append(float(cols[6]))
                                    D33_all.append(float(cols[8]))

                # Process database if entries exist
                if cmul_all:
                    # Convert to numpy arrays
                    cmul_all = np.array(cmul_all)
                    efield_all = np.array(efield_all)
                    Ntz_all = np.array(Ntz_all)
                    Nxi_all = np.array(Nxi_all)
                    D11_all = np.array(D11_all)
                    D31_all = np.array(D31_all)
                    D33_all = np.array(D33_all)

                    # Sort by (efield, cmul) so that e=0 and small cmul come first
                    sort_idx = np.lexsort((efield_all,-cmul_all))
                    cmul_sorted = cmul_all[sort_idx]
                    efield_sorted = efield_all[sort_idx]
                    Nxi_sorted = Nxi_all[sort_idx]
                    Ntz_sorted = Ntz_all[sort_idx]
                    D11_sorted = D11_all[sort_idx]
                    D31_sorted = D31_all[sort_idx]
                    D33_sorted = D33_all[sort_idx]

                    # --- Replace D33 by its efield ~ 0 value for each cmul, keeping array size ---
                    D33_new = D33_sorted.copy()

                    for cm in np.unique(cmul_sorted):
                        mask = (cmul_sorted == cm)
                        idx_min = np.argmin(np.abs(efield_sorted[mask]))
                        D33_at_0 = D33_sorted[mask][idx_min]
                        D33_new[mask] = D33_at_0

                    D33_sorted = D33_new

                    # Compute g11_0 using efield=0 branch and small cmul
                    mask = (efield_sorted == 0.0) & (cmul_sorted >= 1e-5) & (cmul_sorted <= 1e-3)
                    cmul_sel = cmul_sorted[mask]
                    D11_sel = D11_sorted[mask]
                    g11_0 = np.mean(cmul_sel * D11_sel)

                    # Effective ripple
                    eps_eff=(5*g11_0*R_val*R_val*B_val*B_val)**(2./3.)

                    # Range for fit
                    efield_u=0.008*r_val*B_val
                    mask = (cmul_sorted >= 1e-5) & (cmul_sorted <= 1e-3) & \
                           (efield_sorted >= 0) & (efield_sorted <= efield_u)
                    cmul_fit, efield_fit, D11_fit = cmul_sorted[mask], efield_sorted[mask], D11_sorted[mask]

                    # Model for fitting MONKES D11 drift-kernel
                    def D11_model(xdata, g11_er, exp):
                        cmul, efield = xdata
                        term1 = (cmul / g11_0)**exp
                        term2 = ( efield**1.5 / (np.sqrt(cmul) * g11_er ))**exp
                        return (term1 + term2)**(-1/exp)

                    # Fit coefficients
                    p0 = [51E-05, 2.0] # initial guess
                    xdata = (cmul_fit, efield_fit)
                    try:
                        popt, _ = curve_fit(D11_model, xdata, D11_fit, p0=p0, maxfev=10000)
                    except RuntimeError as e:
                        print("Warning: curve_fit not converged",e)
                        popt = p0   # default value
                    g11_er, ex_er = popt 

                    g11_0_scaled    = g11_0 * B_val * B_val
                    eps_eff_scaled  = eps_eff * B_val**(4./3.)
                    g11_er_scaled   = g11_er * np.sqrt(B_val)
                    efield_u_scaled = efield_u / B_val
                    # Write fitted ripple coefficients
                    f.write("c eps_eff g11_ft efield_u g11_er ex_er\n")
                    f.write("cfit {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}\n"
                            .format(eps_eff, g11_0_scaled, efield_u_scaled, g11_er_scaled, ex_er))

                    # Write monoenergetic transport table
                    for cm, ef, Ntz, Nxi, d11, d31, d33 in zip(cmul_sorted, efield_sorted,Ntz_sorted,Nxi_sorted,D11_sorted,D31_sorted,D33_sorted):
                        d11_scaled = -d11 * B_val * B_val
                        d31_scaled = d31 * B_val
                        d33_scaled = -d33
                        ef_scaled = ef / B_val
                        f.write(f"{cm:14.6e} {ef_scaled:14.6e} {d11_scaled:14.6e} "
                                f"{d31_scaled:14.6e} {d33_scaled:14.6e}\n")
                        f.write(" >3{:10d}{:6d}{:14.6E}{:14.6E}{:14.6E}\n".format(int(Ntz),int(Nxi),0.0,0.0,0.0))

                    f.write(" e\n")
                    print(f"Written {len(cmul_sorted)} MONKES entries for surface {surf_idx}")
                else:
                    # No database found for this surface
                    print(f"No MONKES database found for surface {surf_idx}")


