## Overview
Generates cubic turbulent driving fields using the [Mac Low (1999)](https://iopscience.iop.org/article/10.1086/307784) method. Supports solenoidal fraction `f_sol` (1=solenoidal, 0=compressive) and outputs ORION2-ready HDF5 files.

## Dependencies
- FFTW3 (double precision)
- HDF5 C & HL libraries (C API only; no C++ interface)

Homebrew example (Apple Silicon/Intel):
```bash
brew install fftw hdf5
```

## Build
```bash
cmake -S . -B build 
cmake --build build
```

If your paths differ, override `FFTW3_LIB_PATH` / `HDF5_C_LIB_PATH` / `HDF5_HL_LIB_PATH` at configure time.

## Run
```bash
./build/turbdrive
```
Provide when prompted:
1. `n`: grid points per dimension
2. `k`: dimensionless wavenumber (fills shell [k-1, k])
3. `seed`: random seed (`mt19937`, deterministic for repeatability)
4. `mode`: mixing mode (1 or 2; see below)
5. `f_sol`: solenoidal weight (definition differs by mode; see below)

Mode (prompted at runtime):
- Mode 1 — single-field Helmholtz mix (Federrath-style):
  - Build one Gaussian random field in k-space, Helmholtz decomposition into solenoidal and compressive parts, then mix amplitudes with `v = f_sol * v_sol + (1 - f_sol) * v_comp`.
  - Same underlying field: sol has 2 DoF, comp has 1 DoF. Resulting solenoidal energy fraction tends to $\frac{2 f_{sol}^2}{2 f_{sol}^2 + (1 - f_{sol})^2}$. 
  - **Notice:** Here `f_sol` is the $\zeta$ in [Federrath et al. 2010](https://www.aanda.org/articles/aa/full_html/2010/04/aa12437-09/aa12437-09.html) eq. (6). Set `f_sol=0.5` for the so called "natural mixing" (~2/3 solenoidal energy fraction).
- Mode 2 — two-field normalized mix (Similar to Tools/Python/perturbation.py in ORION2 repo, but define f_sol directly the energy fraction):
  - Build two Gaussian random fields, Helmholtz decomposition each; keep only sol from one (seed) and only comp from the other (seed+1), normalize both to rms = 1, then mix in real space with `v = sqrt(f_sol) * v_sol + sqrt(1 - f_sol) * v_comp`.
  - **Notice:** Here the user-specified `f_sol` directly equals the solenoidal energy fraction (components are pre-normalized to equal rms). Set `f_sol=0.6667` directly for the so called "natural mixing".

## Workflow
1) Allocate three k-space arrays (vx, vy, vz) initialized to zero.
2) For each mode in k-space: compute `k_mag`; skip if `k_mag < k-1` or `k_mag > k`, limiting all modes in [k-1, k].
3) Draw amplitudes from a Gaussian around unity `N(1,1)` per component and random phases in [0, 2pi); build complex coefficients.
4) Helmholtz decomposition into parallel (compressive) and perpendicular (solenoidal) parts.
5) Mix per mode choice: See above Mode 1 and Mode 2 descriptions.
6) Inverse FFT (FFTW backward) to real space to obtain pertx, perty, pertz.
7) Save to HDF5 (C API) with datasets `pertx`, `perty`, `pertz` and attributes `n`, `k`, `seed`, `f_sol`, `mode`; filename encodes parameters.

Output file name:
```
zdrv_n{n}_k{k}_seed{seed}_f{f_sol}_mode{mode}.hdf5
```
Datasets: `pertx`, `perty`, `pertz`, shape `(n, n, n)`; attributes: `n`, `k`, `seed`, `f_sol`, `mode`.

## Tools
- Power spectrum: `python Tools/power_spectrum.py /path/to/zdrv_*.hdf5`
  - Loads the perturbation field, computes the 1D power spectrum, and saves `PS_<input_stem>.png` alongside the script.
- Solenoidal/compressive ratio: `python Tools/compute_pert_ratio.py /path/to/zdrv_*.hdf5`
  - FFT-based Helmholtz decomposition; prints solenoidal and compressive energy fractions and their sum check.
- Velocity magnitude contour: `python Tools/pert_contour.py /path/to/zdrv_*.hdf5 [--z ZINDEX]`
  - Contour plot of |v| on a chosen z-slice (default mid-plane); optional `--z` to pick the slice.
- Velocity quiver slice: `python Tools/pert_visual.py /path/to/zdrv_*.hdf5 [--z ZINDEX] [--step N]`
  - Quiver plot of (v_x, v_y) on a z-slice (default mid-plane); `--step` controls downsampling of arrows.

## Status
Still under testing; use cautiously for scientific results.
