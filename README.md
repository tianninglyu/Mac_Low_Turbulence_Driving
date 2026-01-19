## Overview
Generates cubic turbulent driving fields using the Mac Low (1999) method. Supports solenoidal energy fraction `f_sol` (1=solenoidal, 0=compressive) and outputs ORION2-ready HDF5 files.

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
2. `k`: dimensionless wavenumber (fills shell (k-1, k])
3. `seed`: random seed (`mt19937`, deterministic for repeatability)
4. `f_sol`: solenoidal energy fraction (1=solenoidal, 0=compressive)

Output file name:
```
zdrv_n{n}_k{k}_seed{seed}_f{f_sol}.hdf5
```
Datasets: `pertx`, `perty`, `pertz`, shape `(n, n, n)`; attributes: `n`, `k`, `seed`, `f_sol`.

## Tools
- Power spectrum: `python Tools/power_spectrum.py /path/to/zdrv_*.hdf5`
  - Loads the perturbation field, computes the 1D power spectrum, and saves `PS_<input_stem>.png` alongside the script.
- Solenoidal/compressive ratio: `python Tools/compute_pert_ratio.py /path/to/zdrv_*.hdf5`
  - FFT-based Helmholtz projection; prints solenoidal and compressive energy fractions and their sum check.

## Status
Still under testing; use cautiously for scientific results.
