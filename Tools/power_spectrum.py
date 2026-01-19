import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


def calculate_and_plot_power_spectrum(file_path: str) -> None:
    """Read a 3D perturbation field from HDF5 and plot its 1D power spectrum."""
    try:
        with h5py.File(file_path, "r") as f:
            print(f"Loaded file: {file_path}")
            pertx = f["pertx"][:]
            perty = f["perty"][:]
            pertz = f["pertz"][:]
            n = int(f.attrs["n"])
    except (FileNotFoundError, KeyError) as e:
        print(f"Error: cannot read file '{file_path}'. Check path and contents.")
        print(f"Details: {e}")
        return

    pertx_k = np.fft.fftn(pertx, norm="ortho")
    perty_k = np.fft.fftn(perty, norm="ortho")
    pertz_k = np.fft.fftn(pertz, norm="ortho")

    power_k_3d = np.abs(pertx_k) ** 2 + np.abs(perty_k) ** 2 + np.abs(pertz_k) ** 2

    k_freq = np.fft.fftfreq(n) * n
    kx, ky, kz = np.meshgrid(k_freq, k_freq, k_freq, indexing="ij")

    k_magnitude = np.sqrt(kx**2 + ky**2 + kz**2)

    k_bins = np.arange(-0.5, 4.5, 1.0)
    power_spectrum_1d, bin_edges = np.histogram(
        k_magnitude.flatten(), bins=k_bins, weights=power_k_3d.flatten()
    )
    k_values_1d = (bin_edges[:-1] + bin_edges[1:]) / 2

    print("\n k | Power")
    print("-----------")
    for k_val, p in zip(k_values_1d, power_spectrum_1d):
        if k_val % 1 == 0 and k_val >= 0:
            print(f"{int(k_val):>2} | {p:.4e}")

    plt.figure(figsize=(10, 7))
    plt.plot(k_values_1d, power_spectrum_1d, "o-", label="Calculated Spectrum")
    plt.xlabel("Wavenumber, $k$")
    plt.ylabel("Power, $P(k)$")
    plt.title("1D Power Spectrum of 3D Perturbation Field")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

    stem = Path(file_path).stem
    output_filename = f"PS_{stem}.png"
    plt.savefig(output_filename, dpi=150)
    print(f"\nSaved plot to: {output_filename}")
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute 1D power spectrum from HDF5 perturbation field.")
    parser.add_argument("hdf5_file_path", help="Path to input HDF5 file")
    args = parser.parse_args()
    calculate_and_plot_power_spectrum(args.hdf5_file_path)


if __name__ == "__main__":
    main()
