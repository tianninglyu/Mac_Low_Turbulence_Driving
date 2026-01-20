import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt


def plot_quiver(fname: Path, z_index: int | None = None, step: int = 8) -> None:
    with h5py.File(fname, "r") as f:
        pertx = f["pertx"][:]
        perty = f["perty"][:]
        pertz = f["pertz"][:]

    nz = pertx.shape[2]
    zmid = nz // 2 if z_index is None else z_index

    vx = pertx[:, :, zmid]
    vy = perty[:, :, zmid]

    x = np.arange(0, vx.shape[0], step)
    y = np.arange(0, vx.shape[1], step)
    X, Y = np.meshgrid(x, y, indexing="ij")

    U = vx[::step, ::step]
    V = vy[::step, ::step]

    plt.figure(figsize=(6, 6))
    plt.quiver(X, Y, U, V, angles="xy", scale_units="xy")
    plt.xlabel("x index")
    plt.ylabel("y index")
    plt.title(f"Perturbation vector field at z = {zmid}")
    plt.tight_layout()
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Quiver plot of perturbation field slice from HDF5")
    parser.add_argument("file", type=Path, help="path to zdrv_*.hdf5")
    parser.add_argument("--z", type=int, dest="z", default=None, help="z index to plot (default: mid-plane)")
    parser.add_argument("--step", type=int, dest="step", default=8, help="downsample step for quiver arrows")
    args = parser.parse_args()

    plot_quiver(args.file, args.z, args.step)


if __name__ == "__main__":
    main()
