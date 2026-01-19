#!/usr/bin/env python3
"""Compute solenoidal and compressive energy fractions of a turbulence driving field.

Usage:
    python compute_turbulence_ratio.py <path_to_file.hdf5>

Example:
    python compute_turbulence_ratio.py Mac_Low_Turbulence_Driving/zdrv_n512_k2_seed2025.hdf5
"""

import h5py
import numpy as np
import os
import sys

def compute_solenoidal_compressive_ratio(hdf5_path):
    if not os.path.exists(hdf5_path):
        raise FileNotFoundError(f"File not found: {hdf5_path}")

    print(f"Loading file: {hdf5_path}")
    with h5py.File(hdf5_path, "r") as f:
        vx = f["pertx"][:]
        vy = f["perty"][:]
        vz = f["pertz"][:]
        n = f.attrs.get("n", vx.shape[0])
        k_inj = f.attrs.get("k", None)
        seed = f.attrs.get("seed", None)

    print(f"Shape: {vx.shape}, n={n}, k={k_inj}, seed={seed}")

    # --- 3D FFT ---
    vxk = np.fft.fftn(vx)
    vyk = np.fft.fftn(vy)
    vzk = np.fft.fftn(vz)

    # --- Wavenumber grid ---
    kx = np.fft.fftfreq(n) * n
    ky = np.fft.fftfreq(n) * n
    kz = np.fft.fftfreq(n) * n
    kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')

    k2 = kx**2 + ky**2 + kz**2
    k2[0, 0, 0] = 1.0  # avoid division by zero at k=0

    # --- Unit wavevector ---
    kx_hat = kx / np.sqrt(k2)
    ky_hat = ky / np.sqrt(k2)
    kz_hat = kz / np.sqrt(k2)

    # --- Project onto compressive mode: (k·v) k_hat ---
    kv_dot = kx_hat * vxk + ky_hat * vyk + kz_hat * vzk
    vxk_comp = kv_dot * kx_hat
    vyk_comp = kv_dot * ky_hat
    vzk_comp = kv_dot * kz_hat

    # --- Solenoidal mode = total - compressive ---
    vxk_sol = vxk - vxk_comp
    vyk_sol = vyk - vyk_comp
    vzk_sol = vzk - vzk_comp

    # --- Energy ---
    E_comp = np.sum(np.abs(vxk_comp)**2 + np.abs(vyk_comp)**2 + np.abs(vzk_comp)**2)
    E_sol = np.sum(np.abs(vxk_sol)**2 + np.abs(vyk_sol)**2 + np.abs(vzk_sol)**2)
    E_tot = E_comp + E_sol

    f_sol = E_sol / E_tot
    f_comp = E_comp / E_tot

    print("\n--- Decomposition ---")
    print(f"Solenoidal fraction:   {f_sol:.4f}")
    print(f"Compressive fraction:  {f_comp:.4f}")
    print(f"Check: f_sol + f_comp = {f_sol + f_comp:.4f}\n")

    return f_sol, f_comp


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    compute_solenoidal_compressive_ratio(path)

