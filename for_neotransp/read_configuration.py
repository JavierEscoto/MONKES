#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import glob
import os
import subprocess
import sys
import re
from scipy.optimize import curve_fit
from pathlib import Path


# -----------------------------
# Reading magnetic configuration
# -----------------------------
# This function reads either *.bc, boozer.txt, boozmn.nc, ddkes2.data 
# It returns all physical quantities needed to construct MONKES inputs and tables.
def read_configuration(filename: str):

    # Validate input
    if not filename or not isinstance(filename, str):
        raise ValueError("You must provide the file name as a non-empty string.")
    filename = os.path.abspath(os.path.expanduser(os.path.expandvars(filename)))
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    # Check extension
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    
    # -------------------------------
    # 1) Read .bc file (IPP format)
    # ------------------------------
    if ext in (".txt", ".bc"):
        print("Reading data from .bc file")
        with open(filename, "r") as f:
            lines = f.readlines()

        # Extract <beta> from the header
        averagebeta = np.nan
        for line in lines:
            if line.strip().startswith("CC <beta>"):
                try:
                    averagebeta = float(line.strip().split('=')[1])
                except:
                    pass
                break

        # Locate header line containing m0b
        header_index = None
        for line in lines:
            if line.strip().startswith("m0b"):
                header_index = lines.index(line)
                break
        if header_index is None:
            raise ValueError("Cannot find m0b line in boozer.txt")

        # Extract general geometric and magnetic parameters
        values_line = lines[header_index + 1].split()
        m0b, n0b, nsurf, nper = map(int, values_line[:4])
        flux = float(values_line[4])
        minorradius = float(values_line[5])
        majorradius = float(values_line[6])
        minorradiusVMEC = float(values_line[7])
        majorradiusVMEC = float(values_line[8])
        volumeVMEC = float(values_line[9])
        Nperiods = nper

        # Allocate arrays for Fourier components and iota
        s = np.zeros(nsurf)
        iota_arr = np.zeros(nsurf)
        B00 = np.zeros(nsurf)
        B10 = np.zeros(nsurf)
        B01 = np.zeros(nsurf)
        B11 = np.zeros(nsurf)
        R00 = np.zeros(nsurf)
        FSAgpsipsioverB2 = np.zeros(nsurf)
        FSAgpsipsi       = np.zeros(nsurf)
        FSAsqrtgpsipsi   = np.zeros(nsurf)
        FSAinvB2         = np.zeros(nsurf) 
        

        # Loop through surfaces and read Fourier components
        current_line = header_index + 2
        Baxis_phi0 = 0.0
        for surf in range(nsurf):
            # Find dp/ds line
            while "dp/ds" not in lines[current_line]:
                current_line += 1
                if current_line >= len(lines):
                    raise ValueError("Reached end of file while searching for dp/ds")
            current_line += 1
            surf_line = lines[current_line].split()
            s[surf] = float(surf_line[0])
            iota_arr[surf] = float(surf_line[1])
            current_line += 2

            # Read bmn and Rmn Fourier coefficients
            mn_list = []
            while current_line < len(lines):
                cols = lines[current_line].split()
                if len(cols) < 6 or cols[0].lower() == 'm':
                    break
                try:
                    m = int(cols[0])
                    n = int(cols[1])
                    Rmn = float(cols[2])
                    Zmn = float(cols[3])
                    Pmn = -float(cols[4])*2.0*np.pi/Nperiods
                    Bmn = float(cols[5])
                    mn_list.append((m, n, Rmn, Zmn, Pmn, Bmn))
                except ValueError:
                    pass
                current_line += 1

            # Assign Fourier coefficients
            for m, n, Rmn, Zmn, Pmn, Bmn in mn_list:
                if m == 0 and n == 0:
                    B00[surf] = Bmn
                    R00[surf] = Rmn
                elif m == 1 and n == 0:
                    B10[surf] = Bmn
                elif m == 0 and n == 1:
                    B01[surf] = Bmn
                elif m == 1 and n == 1:
                    B11[surf] = Bmn
            # Compute |∇psi|^2/B^2 on a (theta, zeta) grid and its surface-average
            n_theta, n_zeta = 64, 64       # grid resolution (adjust if needed)
            theta = np.linspace(0.0, 2.0*np.pi, n_theta, endpoint=False)
            zeta  = np.linspace(0.0, 2.0*np.pi, n_zeta,  endpoint=False)
            TH, ZE = np.meshgrid(theta, zeta, indexing='ij')
            # Accumulate Fourier sums
            R              = np.zeros_like(TH)
            dR_dtheta      = np.zeros_like(TH)
            dAz_dtheta     = np.zeros_like(TH)
            dZ_dtheta      = np.zeros_like(TH)
            B              = np.zeros_like(TH)
            for m, n, Rmn, Zmn, Pmn, Bmn in mn_list:
                phase = m*TH - Nperiods*n*ZE
                cos_phase = np.cos(phase)
                sin_phase = np.sin(phase)
                R          += Rmn * cos_phase
                dR_dtheta  += (-m) * Rmn * sin_phase
                dAz_dtheta += (-m) * Pmn * cos_phase
                dZ_dtheta  += (+m) * Zmn * cos_phase
                B          += Bmn * cos_phase
            # Compute |∇ψ|^2/B^2
            B2 = B*B
            denominator = np.sum(1.0 / B2)
            nablapsi2_over_B2 = (dR_dtheta**2) + (R**2) * (dAz_dtheta**2) + (dZ_dtheta**2)
            # Surface averages
            FSAgpsipsioverB2[surf] = np.sum(nablapsi2_over_B2 / B2 ) / denominator
            FSAgpsipsi[surf] = np.sum(nablapsi2_over_B2 ) / denominator
            nablapsi_over_B2 = np.sqrt(nablapsi2_over_B2 * B2 ) / B2
            FSAsqrtgpsipsi[surf] = np.sum(nablapsi_over_B2 ) / denominator
            B4 = B2*B2
            FSAinvB2[surf] = np.sum(1.0 / B4)/ denominator

            # Compute on-axis field
            if surf == 0:
                for m, n, Rmn, Zmn, Pmn, Bmn in mn_list:
                    if m == 0:
                        Baxis_phi0 += Bmn

        # Compute radius arrays
        rho = np.sqrt(s)
        rVMEC = rho * minorradiusVMEC

        # Flux normalization
        psi_a = -flux / (2*np.pi)
        
#        write_bc_scaled_copy(src_file=filename, dst_file="boozer_1T.txt", Bscale=1./B00[0])

    # -------------------------------
    # 2) Read .nc format (output of BOOZER_XFORM)
    # ------------------------------
    elif ext in (".nc"):
        print("Reading data from booz_xform file")

        nc = Dataset(filename, "r")

        # From Fortran: nfp_b
        Nperiods = int(nc.variables["nfp_b"][:])

        # Aspect ratio → minor radius
        Aspect_ratio = float(nc.variables["aspect_b"][:])

        # Number of surfaces in VMEC
        nsurf = int(nc.variables["ns_b"][:])-1

#        # Surfaces used in Boozer transform
#        i_b = nc.variables["jlist"][:] - 1  # convert to 0-based index
#        Ns_b = len(i_b)

        # Arrays to be recut (ignore first element)
        iota_arr = nc.variables["iota_b"][1:]
        pres_s = nc.variables["pres_b"][1:]
        beta_s = nc.variables["beta_b"][1:]
        phip_s = nc.variables["phip_b"][1:]
        psi_s = nc.variables["phi_b"][1:]
        B_zeta_s = nc.variables["bvco_b"][1:]
        B_theta_s = nc.variables["buco_b"][1:]

        # Normalized radial coordinate
        s = psi_s / psi_s[-1]

        # Fourier modes
        mnboz_b = int(nc.variables["mnboz_b"][:])
        ixm = nc.variables["ixm_b"][:]
        ixn = nc.variables["ixn_b"][:]
        ixn = ixn / Nperiods  # project to one field period
        bmnc = nc.variables["bmnc_b"][:]  # (mnboz_b, Ns_b)
        rmnc_b = nc.variables["rmnc_b"][:]  # (mnboz_b, Ns_b)
        zmns_b = nc.variables["zmns_b"][:]  # (mnboz_b, Ns_b)
        pmns_b = nc.variables["pmns_b"][:]  # (mnboz_b, Ns_b)

        # Allocate output arrays
        B00 = np.zeros(nsurf)
        B10 = np.zeros(nsurf)
        B01 = np.zeros(nsurf)
        B11 = np.zeros(nsurf)
        R00 = np.zeros(nsurf)
        FSAgpsipsioverB2 = np.zeros(nsurf)
        FSAgpsipsi       = np.zeros(nsurf)
        FSAsqrtgpsipsi   = np.zeros(nsurf)
        FSAinvB2         = np.zeros(nsurf)

        Baxis_phi0=0
        for surf in range(nsurf):
            mn_list = []
            for k in range(mnboz_b):
                m = ixm[k]
                n = ixn[k]
                Rmn = rmnc_b[surf,k]
                Zmn = zmns_b[surf,k] 
                Pmn = pmns_b[surf,k]
                Bmn = bmnc[surf,k]
                mn_list.append((m, n, Rmn, Zmn, Pmn, Bmn))

                if m == 0 and n == 0:
                    B00[surf] = bmnc[surf,k]
                    R00[surf] = rmnc_b[surf,k]
                elif m == 1 and n == 0:
                    B10[surf] = bmnc[surf,k]

                elif m == 0 and n == 1:
                    B01[surf] = bmnc[surf,k]
                elif m == 1 and n == 1:
                    B11[surf] = bmnc[surf,k]
                if surf == 0 and m == 0:
                    Baxis_phi0 += bmnc[surf,k]

            n_theta, n_zeta = 64, 64       # grid resolution (adjust if needed)
            theta = np.linspace(0.0, 2.0*np.pi, n_theta, endpoint=False)
            zeta  = np.linspace(0.0, 2.0*np.pi, n_zeta,  endpoint=False)
            TH, ZE = np.meshgrid(theta, zeta, indexing='ij')
            # Accumulate Fourier sums
            R              = np.zeros_like(TH)
            dR_dtheta      = np.zeros_like(TH)
            dAz_dtheta     = np.zeros_like(TH)  # Σ (?~H~Rm Pmn cos(...)), multiplied by R later
            dZ_dtheta      = np.zeros_like(TH)
            B              = np.zeros_like(TH)
            for m, n, Rmn, Zmn, Pmn, Bmn in mn_list:
                phase = m*TH - Nperiods*n*ZE
                cos_phase = np.cos(phase)
                sin_phase = np.sin(phase)
                R          += Rmn * cos_phase
                dR_dtheta  += (-m) * Rmn * sin_phase
                dAz_dtheta += (-m) * Pmn * cos_phase
                dZ_dtheta  += (+m) * Zmn * cos_phase
                B          += Bmn * cos_phase
            # Compute |?~H~G?~H|^2/B^2
            B2 = B*B
            denominator = np.sum(1.0 / B2)
            nablapsi2_over_B2 = (dR_dtheta**2) + (R**2) * (dAz_dtheta**2) + (dZ_dtheta**2)
            # Surface averages
            FSAgpsipsioverB2[surf] = np.sum(nablapsi2_over_B2 / B2 ) / denominator
            FSAgpsipsi[surf] = np.sum(nablapsi2_over_B2 ) / denominator
            nablapsi_over_B2 = np.sqrt(nablapsi2_over_B2 * B2 ) / B2
            FSAsqrtgpsipsi[surf] = np.sum(nablapsi_over_B2 ) / denominator
            B4 = B2*B2
            FSAinvB2[surf] = np.sum(1.0 / B4)/ denominator

        majorradius = R00[1]
        minorradius = majorradius / Aspect_ratio
        minorradiusVMEC = minorradius
        majorradiusVMEC = majorradius
        volumeVMEC = olumeVMEC = 2 * np.pi**2 * majorradiusVMEC * minorradiusVMEC**2
        averagebeta = np.mean(beta_s)

        rho = np.sqrt(s)
        rVMEC = rho * minorradiusVMEC
        psi_a = psi_s[-1] / (2*np.pi)

        nc.close()
        print("Converted boozmn.nc to Boozer-format quantities.")

#        write_boozmn_scaled_copy(src_file=filename, dst_file="boozmn_1T.nc", Bscale=1./B00[0])

    # -------------------------------
    # 3) Read ddkes2.data format (intermediate output of DKES)
    # -------------------------------
    elif ext in (".data"):
        print("Reading data from ddkes2.data-like file")

        # -----------------------------------
        # Read DKES NAMELIST-style file
        # -----------------------------------
        nzperiod = lalpha = lfout = nrun = 0
        mpolb = ntorb = ibbi = 0
        cmul = efield = psip = chip = 0.0
        btheta = bzeta = 0.0

        # Maximum dimensions for the DKES mode array (Fortran-style)
        M_max = 500
        N_max = 500

        # Allocate borbi(j,i)
        # In Fortran: borbi(-N_max:N_max, -M_max:M_max)
        # In Python we shift to positive indices
        borbi = np.zeros((2*N_max+1, 2*M_max+1))

        # Read file contents
        with open(filename,"r") as f:
            lines = f.readlines()

        # Lowercase version for simple parsing
        text = "\n".join(lines).lower()

        # Helper to extract single-value entries such as "psip = ..."
        def get_value(key, default=None):
            pattern = key.lower() + r"\s*=\s*([eE0-9\+\-\.]+)"
            m = re.search(pattern, text)
            return float(m.group(1)) if m else default

        # Extract fundamental quantities
        nzperiod = int(get_value("nzperiod",1))
        psip     = get_value("psip",0.0)
        chip     = get_value("chip",0.0)
        btheta   = get_value("btheta",0.0)
        bzeta    = get_value("bzeta",0.0)

        # -----------------------------------
        # Read DKES Fourier modes borbi(j,i)
        # -----------------------------------
        mode_expr = r"borbi\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*=\s*([eE0-9\+\-\.]+)"
        Bmnori= None
        B01_val = 0.0
        B10_val = 0.0
        B11_val = 0.0

        for (j,i,val) in re.findall(mode_expr, text):
            j = int(j)
            i = int(i)
            if abs(j) <= N_max and abs(i) <= M_max:
                # Shift indices to positive Python indexing
                borbi[j+N_max, i+M_max] = float(val)
                # Save specific modes
                if j == 0 and i == 0:
                    B00_val = float(val)
                elif j == 1 and i == 0:
                    B10_val = float(val)
                elif j == 0 and i == 1:
                    B01_val = float(val)
                elif j == 1 and i == 1:
                    B11_val = float(val)
                # Maximum value excluding (0,0)
                if not (j == 0 and i == 0):
                    if Bmnori is None or abs(float(val)) > Bmnori:
                        Bmnori = abs(float(val))
        # -----------------------------------
        # Read 0-D parameters from filename
        # -----------------------------------

        factB0 = 1.0
        factR0 = 1.0
        facta0 = 1.0
        factBmn= 1.0
        s0ori = -1.0

        if Path(filename).name != "ddkes2.data":

            B0ori   = B00_val
            match = re.search(r"B0_([-+]?\d*\.?\d+)_", filename)
            if match:
                B0 = float(match.group(1))
                factB0 = B0 / B0ori

            R0ori = bzeta / B0ori
            match = re.search(r"R0_([-+]?\d*\.?\d+)_", filename)
            if match:
                R0 = float(match.group(1))
                factR0 = R0 / R0ori
        
            a0ori = abs(psip  / B0ori)
            match = re.search(r"a0_([-+]?\d*\.?\d+)_", filename)
            if match:
                a0 = float(match.group(1))
                facta0 = a0/ a0ori

            s0ori = 1.0
            match = re.search(r"s0_([-+]?\d*\.?\d+).", filename)
            if match:
                s0ori = float(match.group(1))

            Bmn = B0ori * factB0 * np.sqrt( s0ori ) * a0ori * facta0 / ( R0ori * factR0 )
            factBmn= Bmn / Bmnori

        # -----------------------------------
        # Reconstruct equilibrium-like quantities from DKES input
        # -----------------------------------

        # Field periodicity
        Nperiods = nzperiod

        # Rotational transform from DKES parameters
        psi_p = psip * factB0 * facta0  
        chi_p = chip * factB0 * facta0
        iota_val = chi_p / psi_p 

        # Contravariant components (B · ∇theta, B · ∇zeta)
        B_theta = btheta * factB0 * factR0
        B_zeta  = bzeta  * factB0 * factR0

        # Main Fourier modes
        B00_val = B00_val * factB0
        B10_val = B10_val * factBmn
        B01_val = B01_val * factBmn
        B11_val = B11_val * factBmn 

        # -----------------------------------
        # Construct minimal arrays (nsurf = 10)
        # -----------------------------------
        nsurf = 9
        s = np.linspace(0.1,1,nsurf)
        rho = np.sqrt(s)

        iota_arr = np.full(nsurf, iota_val)
        B00 = np.full(nsurf, B00_val)
        B10 = np.full(nsurf, B10_val)
        if s0ori>=0:
            scaling=rho/np.sqrt(s0ori)
        else:
            scaling=rho/rho
        B01 = B01_val * scaling
        B11 = B11_val * scaling
        Baxis_phi0 = B00_val+B10_val

        FSAgpsipsioverB2=np.full(nsurf,np.nan)
        FSAgpsipsi      =np.full(nsurf,np.nan)
        FSAsqrtgpsipsi  =np.full(nsurf,np.nan)
        FSAinvB2        =np.full(nsurf,np.nan)

        # -----------------------------------
        # Large-aspect-ratio approximations (as in the Fortran code)
        # -----------------------------------
        majorradius = abs(B_zeta / B00_val)
        minorradius = abs(psi_p  / B00_val)
        majorradiusVMEC = majorradius
        minorradiusVMEC = minorradius
        volumeVMEC = 2 * np.pi**2 * majorradiusVMEC * minorradiusVMEC**2
        rVMEC = rho * minorradiusVMEC
        R00 = np.full(nsurf,majorradius)

        # Normalization of flux
        psi_a = psi_p * minorradius / 2

        # No beta information in DKES → set to zero
        averagebeta = 0.0

        print("Converted ddkes2.data to Boozer-format quantities.")

        for i in range(1, 8):
            s_i = ((i - 0.5)**2) / (7.0**2)
            if i == 1:
                s_i *= 3
            folder_s = f"s_{s_i:.3f}"
            if os.path.exists(folder_s):
                break
            os.makedirs(folder_s, exist_ok=True)
            #dst_file = os.path.join(folder_s,f"ddkes2_s_{s_i:.3f}.data")
            dst_file = os.path.join(folder_s,f"ddkes2.data")
            facts0=1
            if s0ori > 0:
                facts0 = s_i / s0ori
            print("scaling",filename,factB0,factBmn,facta0,facts0)
            write_ddkes_scaled_copy(filename,dst_file,factB0,factBmn,factR0,facta0,facts0)
#			with open(dst_file, "w") as f:
#				f.write("&datain\n")
#				f.write(f"nzperiod= 5,\n")
#				f.write(f"lalpha= 20, lfout= 20, nrun= 1,\n")
#				f.write(f"cmul= 0.00001,\n")
#				f.write(f"efield= 0.00000,\n")
#				f.write(f"mpolb= 32, ntorb= 64, ibbi= 1,\n")
#				f.write(f"chip= {chip:.6f}, psip= {psip:.6f}, btheta= {btheta:.6f}, bzeta= {bzeta:.6f},\n")
#				for j_shifted in range(borbi.shape[0]):
#					for i_shifted in range(borbi.shape[1]):
#						val = borbi[j_shifted, i_shifted]
#						if val != 0.0:
#							j = j_shifted - N_max
#							i = i_shifted - M_max
#							f.write(f"borbi({j},{i})= {val:.8f},\n")
#				f.write("/\n")
    

    # -------------------------------
    # 4) If neither file exists → fallback to wout.nc
    # -------------------------------
#    elif ext in (".data"):
#        print("boozer.txt and boozmn.nc not found, reading wout.nc")
#        wout = Dataset(wout_path, "r")
#
#        phipf = wout.variables['phipf'][:]
#        psi_a = phipf[0] / (2*np.pi)
#        averagebeta = wout.variables['betatotal'][:]
#        volumeVMEC = wout.variables['volume_p'][:]
#        minorradiusVMEC = wout.variables['Aminor_p'][:]
#        majorradiusVMEC = wout.variables['Rmajor_p'][:]
#        Nperiods = wout.variables['nfp'][:]
#        nsurf = wout.variables['ns'][:]
#        iota_arr = wout.variables['iotaf'][:]
#        bmnc = wout.variables['bmnc'][:]
#
#        s = np.linspace(0,1,nsurf)
#        rho = np.sqrt(s)
#        rVMEC = rho * minorradiusVMEC
#        B00 = np.max(bmnc, axis=1)
#
#        B10 = np.full(nsurf, np.nan)
#        B01 = np.full(nsurf, np.nan)
#        B11 = np.full(nsurf, np.nan)
#        R00 = np.full(nsurf, np.nan)
#        Baxis_phi0 = np.nan
#
#        wout.close()

    else:
        raise FileNotFoundError("No appripriate equilibrium file provided")

    # -------------------------------
    # Return dictionary
    # -------------------------------
    return dict(
        averagebeta=averagebeta,
        s=s, rho=rho, rVMEC=rVMEC, iota_arr=iota_arr,
        B00=B00, B10=B10, B01=B01, B11=B11, R00=R00,
        psi_a=psi_a, volumeVMEC=volumeVMEC,
        minorradiusVMEC=minorradiusVMEC, majorradiusVMEC=majorradiusVMEC,
        Nperiods=Nperiods, Baxis_phi0=Baxis_phi0,
        FSAgpsipsioverB2=FSAgpsipsioverB2, FSAgpsipsi=FSAgpsipsi,
        FSAsqrtgpsipsi=FSAsqrtgpsipsi, FSAinvB2=FSAinvB2
    )



# ------------------------------------------------------------
# Function: write_boozmn_scaled_copy
# ------------------------------------------------------------
def write_boozmn_scaled_copy(src_file="boozmn.nc", dst_file="boozmn_1T.nc", Bscale=1.0):
    """
    Copy all variables from src_file to dst_file, scaling selected ones by Bscale.
    Creates dst_file only if it does not exist.
    """
    if os.path.exists(dst_file):
        print(f"{dst_file} already exists. Nothing done.")
        return

    # Variables to scale
    scale_map = {
        "pres_b": Bscale**2,
        "phip_b": Bscale,
        "phi_b": Bscale,
        "bvco_b": Bscale,
        "buco_b": Bscale,
        "bmnc_b": Bscale,
        "gmn_b": Bscale**-2,
    }

    with Dataset(src_file, "r") as src, Dataset(dst_file, "w") as dst:
        # Copy dimensions
        for name, dim in src.dimensions.items():
            dst.createDimension(name, (len(dim) if not dim.isunlimited() else None))

        # Copy global attributes
        for attr in src.ncattrs():
            dst.setncattr(attr, getattr(src, attr))

        # Copy variables
        for name, var in src.variables.items():
            out_var = dst.createVariable(name, var.dtype, var.dimensions)
            out_var.setncatts({a: var.getncattr(a) for a in var.ncattrs()})

            data = var[...]
            if name in scale_map:
                data = data * scale_map[name]

            out_var[...] = data

    print(f"Created scaled file: {dst_file}")



# ------------------------------------------------------------
# Function: write_ddkes_scaled_copy
# ------------------------------------------------------------
def write_ddkes_scaled_copy(src_file, dst_file, factB0, factBmn, factR0, facta0, facts0):
    """
    Reads ddkes2.data, scales chip, psip, btheta, bzeta and borbi(j,i) by Bscale,
    keeps other variables unchanged, and writes to dst_file.
    """
    with open(src_file, "r") as f:
        text = f.read().lower()

    # Helper to extract values
    def get_int(key, default=0):
        m = re.search(rf"{key}\s*=\s*(\d+)", text)
        return int(m.group(1)) if m else default

    def get_float(key, default=0.0):
        m = re.search(rf"{key}\s*=\s*([eEdD0-9\+\-\.]+)", text)
        return float(m.group(1).replace("d", "e").replace("D", "e")) if m else default

    def fmt(val):  # Fortran-style scientific notation
        return f"{val:.6E}"

    # Unchanged variables
    nzperiod = get_int("nzperiod", 5)
    lalpha   = get_int("lalpha", 20)
    lfout    = get_int("lfout", 20)
    nrun     = get_int("nrun", 1)
    cmul     = get_float("cmul", 3.0)
    efield   = get_float("efield", 0.0)
    mpolb    = get_int("mpolb", 7)
    ntorb    = get_int("ntorb", 31)
    ibbi     = get_int("ibbi", 1)

    # Scaled variables
    chip   = get_float("chip", 0.0) * factB0 * facta0 * np.sqrt(facts0)
    psip   = get_float("psip", 0.0) * factB0 * facta0 * np.sqrt(facts0)
    btheta = get_float("btheta", 0.0) * factB0 * factR0
    bzeta  = get_float("bzeta", 0.0)  * factB0 * factR0

    # Read borbi entries and scale
    borbi_entries = []
    for j, i, val in re.findall(r"borbi\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*=\s*([eEdD0-9\+\-\.]+)", text):
        if int(i) == 0 and int(j) == 0:
            Bscale=factB0
        else:
            Bscale=factBmn
        val_scaled = float(val.replace("D", "E").replace("d", "e")) * Bscale
        if int(i) == 0  and int(j) == 0:
            val_scaled = val_scaled
        elif facts0 < 1:
            val_scaled *= facts0**(float(i)/2)
        else:
            val_scaled *= np.sqrt(facts0)
        borbi_entries.append((int(j), int(i), val_scaled))

    # Write output
    with open(dst_file, "w") as f:
        f.write("&datain\n")
        f.write(f" nzperiod= {nzperiod}, lalpha= {lalpha}, lfout={lfout}, nrun = {nrun},\n")
        f.write(f" cmul = {cmul:.4f}, efield = {efield:.4f},\n")
        f.write(f" mpolb = {mpolb}, ntorb = {ntorb}, ibbi = {ibbi},\n")
        f.write(f" chip = {fmt(chip)}, psip = {fmt(psip)}, btheta = {fmt(btheta)}, bzeta = {fmt(bzeta)},\n")
        for j, i, val in borbi_entries:
            f.write(f" borbi({j},{i})= {fmt(val)},\n")
        f.write("/\n")

    print(f"Created scaled file: {dst_file}")
    return dst_file



# ------------------------------------------------------------
# Function: write_bc_scaled_copy
# ------------------------------------------------------------

def write_bc_scaled_copy(src_file, dst_file, Bscale):

    with open(src_file, "r") as f:
        lines = f.readlines()

    # Copy original lines and create a dictionary for replacements by index
    original_lines = lines[:]
    replacements = {}

    # Locate header line containing m0b
    header_index = None
    for line in lines:
        if line.strip().startswith("m0b"):
            header_index = lines.index(line)
            break

    # Extract general geometric and magnetic parameters
    values_line = lines[header_index + 1].split()
    m0b, n0b, nsurf, nper = map(int, values_line[:4])
    flux = float(values_line[4])
    minorradius = float(values_line[5])
    majorradius = float(values_line[6])
    minorradiusVMEC = float(values_line[7])
    majorradiusVMEC = float(values_line[8])
    volumeVMEC = float(values_line[9])

    # Rewrite the general parameters line with scaled flux
    try:
        _vals = lines[header_index + 1].split()
        _vals[4] = f"{flux * Bscale:.12g}"
        replacements[header_index + 1] = " ".join(_vals) + ("\n" if lines[header_index + 1].endswith("\n") else "")
    except Exception:
        pass

    # Loop through surfaces and read Fourier components
    current_line = header_index + 2
    for surf in range(nsurf):
        # Find the line containing "dp/ds"
        while "dp/ds" not in lines[current_line]:
            current_line += 1
        current_line += 1
        surf_line = lines[current_line].split()
        s = float(surf_line[0])
        iota_arr = float(surf_line[1])

        curr_pol = float(surf_line[2])
        curr_tor = float(surf_line[3])
        pprime = float(surf_line[4])

        current_line += 2

        # Rewrite the surface line with scaled curr_pol, curr_tor, and pprime
        try:
            _surf_vals = lines[current_line - 2].split()  # line with: s iota curr_pol curr_tor pprime
            if len(_surf_vals) >= 5:
                _surf_vals[2] = f"{curr_pol * Bscale:.12g}"
                _surf_vals[3] = f"{curr_tor * Bscale:.12g}"
                _surf_vals[4] = f"{pprime * Bscale * Bscale:.4E}"
                replacements[current_line - 2] = " ".join(_surf_vals) + ("\n" if lines[current_line - 2].endswith("\n") else "")
        except Exception:
            pass

        # Read bmn and Rmn Fourier coefficients
        while current_line < len(lines):
            cols = lines[current_line].split()
            if len(cols) < 6 or cols[0].lower() == 'm':
                break
            try:
                m = int(cols[0])
                n = int(cols[1])
                Rmn = float(cols[2])
                Zmn = float(cols[3])
                pmn = float(cols[4])
                Bmn = float(cols[5])
            except ValueError:
                pass

            try:
                newline = "\n" if lines[current_line].endswith("\n") else ""
                line_out = (
                    f"{m:d} {n:d} "
                    f"{Rmn: .8E} {Zmn: .8E} {pmn: .8E} {Bmn * Bscale: .8E}"
                    f"{newline}"
                )
                replacements[current_line] = line_out
            except Exception:
                pass

            current_line += 1

    # Write all lines to dst_file, applying replacements where needed
    try:
        with open(dst_file, "w") as out:
            for i, line in enumerate(original_lines):
                # ADDED: append extra text to the first line
                if i == n0:
                    line0 = replacements.get(0, line)
                    newline = "\n" if line0.endswith("\n") else ""
                    base = line0[:-1] if line0.endswith("\n") else line0
                    line0_mod = base + "         CC SCALED TO B00=1T AT THE INNERMOST FLUX-SURFACE" + newline
                    out.write(line0_mod)
                else:
                    out.write(replacements.get(i, line))
        print(f"Scaled .bc copy written to: {dst_file}")
    except Exception as e:
        print(f"Failed writing scaled copy: {e}")

