#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import glob
import os
import subprocess
import sys
from scipy.optimize import curve_fit

# -----------------------------
# Write monkes.dk header
# -----------------------------
# This function writes an initial block of equilibrium data to monkes.dk,
# which is later appended with MONKES monoenergetic tables.
def write_header(params, filename="monkes.dk"):
    s = params['s']
    rho = params['rho']
    rVMEC = params['rVMEC']
    iota_arr = params['iota_arr']
    B00 = params['B00']
    B10 = params['B10']
    B01 = params['B01']
    B11 = params['B11']
    psi_a = params['psi_a']
    averagebeta = params['averagebeta']
    minorradiusVMEC = params['minorradiusVMEC']
    majorradiusVMEC = params['majorradiusVMEC']
    volumeVMEC = params['volumeVMEC']
    Nperiods = params['Nperiods']
    Baxis_phi0 = params['Baxis_phi0']
    FSAgpsipsioverB2 = params['FSAgpsipsioverB2']
    FSAgpsipsi = params['FSAgpsipsi']
    FSAsqrtgpsipsi = params['FSAsqrtgpsipsi']
    FSAinvB2 = params['FSAinvB2']
    nsurf = len(s)

    # Volume and its derivatives for each surface
    dVds = np.full(nsurf, volumeVMEC)
    dVdr = dVds * (2*np.sqrt(s)/minorradiusVMEC)

    with open(filename, "w", encoding="utf-8") as f:
        # Write general equilibrium metadata
        f.write(f"cc   <beta> = {averagebeta}\n")
        f.write("cc begin data from bc file\n")
        f.write(f"cc {'psi_a':14} {'averagebeta':14} {'minorradius':14} {'minorradiusVMEC':14} "
                f"{'majorradiusVMEC':14} {'volumeVMEC':14} {'Nperiods':14} {'Baxis_phi0':14}\n")
        f.write(f"cc {psi_a:14.6e} {averagebeta:14.6e} {minorradiusVMEC:14.6e} {minorradiusVMEC:14.6e} "
                f"{majorradiusVMEC:14.6e} {volumeVMEC:14.6e} {Nperiods:14} {Baxis_phi0:14.6e}\n")

        # Write surface-dependent quantities
        f.write("cc s          rho       r         rVMEC     iota      dVds          dVdr          "
                "B00           B10overB00    B01overB00    B11overB00    FSAgpsipsi   FSAgpsipsioverB2 "
                "FSAsqrtgpsipsi FSAinvB2\n")

        for i in range(nsurf):
            # Normalize Fourier modes by B00
            B10_over_B00 = B10[i]/B00[i] if B00[i]!=0 else np.nan
            B01_over_B00 = B01[i]/B00[i] if B00[i]!=0 else np.nan
            B11_over_B00 = B11[i]/B00[i] if B00[i]!=0 else np.nan

            # Write line for each surface
            f.write(f"cc {s[i]:12.6e} {rho[i]:12.6e} {rVMEC[i]:12.6e} {rVMEC[i]:12.6e} "
                    f"{iota_arr[i]:12.6e} {dVds[i]:12.6e} {dVdr[i]:12.6e} "
                    f"{B00[i]:12.6e} {B10_over_B00:12.6e} {B01_over_B00:12.6e} {B11_over_B00:12.6e} "
                    f"{FSAgpsipsi[i]:12.6e} {FSAgpsipsioverB2[i]:12.6e} {FSAsqrtgpsipsi[i]:12.6e} {FSAinvB2[i]:12.6e}\n")

        f.write("cc end data from bc file\n")
        f.write("cc gradient definition: dPsi/dr = 2 psi_a r/a^2\n")
    print(f"{filename} header written.")


