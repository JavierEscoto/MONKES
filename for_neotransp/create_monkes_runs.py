#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import glob
import os
import subprocess
import sys
from scipy.optimize import curve_fit

# -----------------------------
# Loop surfaces and cmul, create folders and files
# -----------------------------
# This function prepares all folders and input files needed for MONKES runs.
# It also automatically submits jobs via sbatch.
def create_monkes_runs(params,filename: str):

    _, ext = os.path.splitext(filename)
    ext = ext.lower()

    s_arr = params['s']
    iota_arr = params['iota_arr']
    minorradiusVMEC = params['minorradiusVMEC']
    majorradiusVMEC = params['majorradiusVMEC']
    Nperiods = params['Nperiods']
    
    cmul = np.array([
        1.0e+1, 3.0e+0, 1.0e+0, 3.0e-1, 1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3,
        1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6, 3.0e-7
    ])

    # Electric field grid for monoenergetic transport scans
    E_r_values = np.array([
        0.0,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,
        2e-3,5e-3,0.01,0.02,0.03,0.05,0.1,0.2,
        0.3,0.5,0.7,0.8,1,1.2,1.5,2,3,5
    ])
    
    # Loop through the 7 surface definitions used by MONKES
    for i in range(1, 8):
        s_i = ((i - 0.5)**2) / (7.0**2)
        if i == 1:
            s_i *= 3
        folder_s = f"s_{s_i:.3f}"

        # Interpolate rotational transform
        i0 = np.interp(s_i, s_arr, iota_arr)
        r_val = minorradiusVMEC * np.sqrt(s_i)

        # Compute E_r for this radius
        E_r = E_r_values * abs(i0 * r_val / majorradiusVMEC)

        # Loop on collisionality
        for nu_val in cmul:
            folder_nu = os.path.join(folder_s, f"cmul_{nu_val:.1e}")
            if os.path.exists(folder_nu):
                break
            os.makedirs(folder_nu, exist_ok=True)

            if ext in (".txt", ".bc"):
                link_name="boozer.txt"
                original_file="../../" + filename
            elif ext in (".nc"):
                link_name="boozmn.nc"
                original_file="../../" + filename
            elif ext in (".data"):
                link_name="ddkes2.data"
                original_file=f"../ddkes2.data"
                #original_file=f"../ddkes2_s_{s_i:.3f}.data"

            link_path = os.path.join(folder_nu, link_name)
            os.symlink(original_file, link_path)

            # Write monkes_input.surface file
            with open(os.path.join(folder_nu, "monkes_input.surface"), "w") as f:
                f.write("&surface\n")
                f.write(f"s={s_i:.6e}\n")
                f.write("/\n")

            # Select resolution depending on collisionality
            if nu_val > 1e-4:
                N_theta, N_zeta, N_xi = 10, 40, 80
            else:
                N_theta, N_zeta, N_xi = 20, 80, 180
            N_theta, N_zeta, N_xi = 120, 120, 200

            # Write monkes_input.parameters
            with open(os.path.join(folder_nu, "monkes_input.parameters"), "w") as f:
                f.write("&parameters\n")
                f.write(f"N_theta = {N_theta}\n")
                f.write(f"N_zeta = {N_zeta}\n")
                f.write(f"N_xi = {N_xi}\n")
                f.write(f"nu = {nu_val:.6e}\n")
                f.write("E_r = " + "  ".join(f"{v:.6e}" for v in E_r) + "\n")
                f.write("/\n")

            # Submit MONKES run to Slurm cluster
            old_dir = os.getcwd()
            os.chdir(folder_nu)
            os.system("sbatch /u/jove/MONKES/for_neotransp/monkes.sbs")
#            os.system("sbatch /u/jove/MONKES/for_neotransp/monkes_long.sbs")
            os.chdir(old_dir)

        if i == 7:
            print("MONKES input folders and scripts created. MONKES is running")
            return
            
    print("MONKES input folders and scripts created.")

