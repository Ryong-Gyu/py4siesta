import glob
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


class BandStructureData:
    def __init__(self, kpath, energies, nbands, nspin, special_k, labels, fermi_level, bandgap, vbm):
        self.kpath = kpath
        self.energies = energies
        self.nbands = nbands
        self.nspin = nspin
        self.special_k = special_k
        self.labels = labels
        self.fermi_level = fermi_level
        self.bandgap = bandgap
        self.vbm = vbm


def _find_bands_file(bands_path=None):
    if bands_path:
        path = Path(bands_path).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Band file does not exist: {path}")
        return path

    matches = sorted(glob.glob("*.bands"))
    if not matches:
        raise FileNotFoundError("No *.bands file found in the current directory.")
    return Path(matches[0])


def _find_pdos_file(pdos_path=None):
    if pdos_path:
        path = Path(pdos_path).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"PDOS file does not exist: {path}")
        return path

    matches = sorted(glob.glob("*.PDOS"))
    if not matches:
        raise FileNotFoundError("No *.PDOS file found in the current directory.")
    return Path(matches[0])


def _clean_k_label(label):
    cleaned = label.strip().strip("'\"")
    if cleaned in {"\\Gamma", "Gamma", "G"}:
        return r"$\Gamma$"
    return cleaned


def read_band_structure(bands_path=None):
    path = _find_bands_file(bands_path)
    lines = path.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError(f"{path} is too short to be a SIESTA .bands file.")

    fermi_level = float(lines[0].split()[0])
    nbands, nspin, nkpoints = [int(value) for value in lines[3].split()[:3]]
    total_bands = nbands * nspin
    lines_per_kpoint = int(math.ceil(float(total_bands) / 10.0))

    kpath = np.zeros(nkpoints, dtype=float)
    energies = np.zeros((total_bands, nkpoints), dtype=float)
    line_index = 4

    for ikpoint in range(nkpoints):
        band_index = 0
        for segment_index in range(lines_per_kpoint):
            words = lines[line_index].split()
            line_index += 1
            if segment_index == 0:
                kpath[ikpoint] = float(words[0])
                values = words[1:]
            else:
                values = words

            for value in values:
                if band_index >= total_bands:
                    break
                energies[band_index, ikpoint] = float(value)
                band_index += 1

    nspecial = int(lines[line_index].split()[0])
    line_index += 1

    special_k = []
    labels = []
    for line in lines[line_index:line_index + nspecial]:
        words = line.split()
        if len(words) < 2:
            continue
        special_k.append(float(words[0]))
        labels.append(_clean_k_label(words[1]))

    below_fermi = energies[energies <= fermi_level]
    above_fermi = energies[energies > fermi_level]
    vbm = float(np.max(below_fermi)) if below_fermi.size else fermi_level
    cbm = float(np.min(above_fermi)) if above_fermi.size else fermi_level
    bandgap = max(0.0, cbm - vbm)

    return BandStructureData(
        kpath=kpath,
        energies=energies,
        nbands=nbands,
        nspin=nspin,
        special_k=np.array(special_k, dtype=float),
        labels=labels,
        fermi_level=fermi_level,
        bandgap=bandgap,
        vbm=vbm,
    )


def plot_band_structure(bands_path=None, emin=-2.0, emax=4.0, output_path="band.png"):
    data = read_band_structure(bands_path)
    energy_reference = data.vbm
    shifted_energies = data.energies - energy_reference

    fig, ax = plt.subplots(figsize=(3.0, 4.0))

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    ax.tick_params(axis="x", direction="in", length=6)
    ax.tick_params(axis="y", direction="in", length=6)

    colors = ["k", "r", "tab:blue", "tab:green"]
    for spin_index in range(data.nspin):
        start = spin_index * data.nbands
        end = start + data.nbands
        color = colors[spin_index % len(colors)]
        for band_index in range(start, end):
            ax.plot(data.kpath, shifted_energies[band_index], linewidth=2, color=color)

    ax.set_xlim(float(np.min(data.kpath)), float(np.max(data.kpath)))
    ax.set_ylim(float(emin), float(emax))
    ax.set_xticks(data.special_k)
    ax.set_xticklabels(data.labels, fontsize=16)
    ax.set_ylabel(r"E-$E_V$", fontsize=16)
    ax.tick_params(axis="y", labelsize=16)

    for kpoint in data.special_k:
        ax.axvline(x=kpoint, color="k", linestyle="--", linewidth=2.0)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, transparent=True)
    plt.close(fig)

    np.savetxt("specialk.csv", data.special_k, delimiter=",")
    np.savetxt("kpath.csv", data.kpath, delimiter=",")
    np.savetxt("band.csv", shifted_energies.T, delimiter=",")

    return {
        "figure": Path(output_path),
        "special_k": Path("specialk.csv"),
        "kpath": Path("kpath.csv"),
        "bands": Path("band.csv"),
        "fermi_level": data.fermi_level,
        "vbm": data.vbm,
        "bandgap": data.bandgap,
    }


def _group_atoms_by_z(atoms, tolerance=0.01):
    z_coords = []
    indices = []

    for atom in atoms:
        z_coord = atom[2]
        for group_index, existing_z in enumerate(z_coords):
            if abs(z_coord - existing_z) < tolerance:
                indices[group_index].append(atom.get_serial())
                break
        else:
            z_coords.append(z_coord)
            indices.append([atom.get_serial()])

    return z_coords, indices


def _read_fermi_level(eig_path):
    with eig_path.open() as eig_file:
        first_line = eig_file.readline().split()
    if not first_line:
        raise ValueError(f"{eig_path} does not contain a Fermi level.")
    return float(first_line[0])


def plot_pldos(
    pdos_path=None,
    emin=-4.0,
    emax=2.0,
    zmin=None,
    zmax=None,
    broad=0.02,
    npoints=1001,
    output_path="pldos.png",
):
    from NanoCore import io, s2

    path = _find_pdos_file(pdos_path).resolve()
    work_dir = path.parent
    label = path.name[:-5] if path.name.endswith(".PDOS") else path.stem
    xyz_path = work_dir / f"{label}.xyz"
    eig_path = work_dir / f"{label}.EIG"

    if not xyz_path.is_file():
        raise FileNotFoundError(f"Required XYZ file does not exist: {xyz_path}")
    if not eig_path.is_file():
        raise FileNotFoundError(f"Required EIG file does not exist: {eig_path}")

    previous_dir = Path.cwd()
    try:
        os.chdir(work_dir)
        atoms = io.read_xyz(xyz_path.name)._atoms
        fermi_level = _read_fermi_level(Path(eig_path.name))
        z_coords, indices = _group_atoms_by_z(atoms)

        projected_dos = []
        energy = None
        for atom_indices in indices:
            energy_values, dos_up, unused_dos_down = s2.get_pdos(
                None,
                emin,
                emax,
                by_atom=1,
                atom_index=atom_indices,
                broad=broad,
                npoints=npoints,
                label=label,
            )
            energy = np.array(energy_values, dtype=float)
            projected_dos.append(np.array(dos_up, dtype=float))

        if energy is None or not projected_dos:
            raise ValueError(f"No PDOS data could be read from {path}.")

        z_values = np.array(z_coords, dtype=float)
        dos_grid = np.abs(np.array(projected_dos, dtype=float)).T
        log_dos = np.log10(dos_grid + 10 ** -5)
        x_grid, y_grid = np.meshgrid(z_values, energy - fermi_level)

        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(1.5)
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)
        ax.tick_params(axis="x", direction="in", length=12)
        ax.tick_params(axis="y", direction="in", length=12)

        z_floor = -4.0
        z_ceiling = 2.0
        levels = np.linspace(z_floor, z_ceiling, 100)
        contour = ax.contourf(
            x_grid,
            y_grid,
            log_dos,
            levels,
            cmap=plt.get_cmap("viridis"),
            extend="both",
        )
        colorbar = fig.colorbar(contour, ticks=np.linspace(z_floor, z_ceiling, 3), extend="both")
        colorbar.ax.tick_params(labelsize=14)

        ax.set_ylim(float(emin), float(emax))
        ax.set_xlim(
            float(np.min(z_values) if zmin is None else zmin),
            float(np.max(z_values) if zmax is None else zmax),
        )
        ax.tick_params(labelbottom=False, labelleft=False)

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, transparent=True)
        plt.close(fig)

        np.savetxt("pldos_z.csv", z_values, delimiter=",")
        np.savetxt("pldos_energy.csv", energy - fermi_level, delimiter=",")
        np.savetxt("pldos.csv", log_dos, delimiter=",")

        return {
            "figure": work_dir / output_path,
            "z": work_dir / "pldos_z.csv",
            "energy": work_dir / "pldos_energy.csv",
            "pldos": work_dir / "pldos.csv",
            "fermi_level": fermi_level,
        }
    finally:
        os.chdir(previous_dir)
