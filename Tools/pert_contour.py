import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt


def plot_slice(fname: Path, z_index: int | None = None) -> None:
    with h5py.File(fname, "r") as f:
        vx = f["pertx"][:]
        vy = f["perty"][:]
        vz = f["pertz"][:]

    zmid = vx.shape[2] // 2 if z_index is None else z_index

    speed = np.sqrt(
        vx[:, :, zmid] ** 2 +
        vy[:, :, zmid] ** 2 +
        vz[:, :, zmid] ** 2
    )

    plt.figure(figsize=(6, 5))
    cs = plt.contourf(
        speed.T,
        levels=40,
        cmap="viridis"
    )
    plt.colorbar(cs, label=r"$|\\mathbf{v}|$")
    plt.xlabel("x index")
    plt.ylabel("y index")
    plt.title(f"Velocity magnitude at z = {zmid}")
    plt.tight_layout()
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Contour plot of velocity magnitude slice from HDF5 perturbation field")
    parser.add_argument("file", type=Path, help="path to zdrv_*.hdf5")
    parser.add_argument("--z", type=int, dest="z", default=None, help="z index to plot (default: mid-plane)")
    args = parser.parse_args()

    plot_slice(args.file, args.z)


if __name__ == "__main__":
    main()
