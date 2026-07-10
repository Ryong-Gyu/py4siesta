import glob
import math
import os
import re
import subprocess
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
    ax.set_ylabel(r"$E-E_V$", fontsize=16)
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


def _find_matching_file(label, suffix, work_dir):
    path = work_dir / f"{label}{suffix}"
    if path.is_file():
        return path

    matches = sorted(work_dir.glob(f"*{suffix}"))
    if not matches:
        raise FileNotFoundError(f"No *{suffix} file found in {work_dir}.")
    return matches[0]


def _find_optional_matching_file(label, suffix, work_dir):
    try:
        return _find_matching_file(label, suffix, work_dir)
    except FileNotFoundError:
        return None


def _read_eig_levels(eig_path):
    lines = eig_path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{eig_path} is too short to be a SIESTA .EIG file.")

    fermi_level = float(lines[0].split()[0])
    neig, nspin, nkpoints = [int(value) for value in lines[1].split()[:3]]
    total_eigenvalues = neig * nspin
    lines_per_kpoint = int(math.ceil(float(total_eigenvalues) / 10.0))

    energies = np.zeros((nkpoints, total_eigenvalues), dtype=float)
    line_index = 2
    for ikpoint in range(nkpoints):
        eigenvalue_index = 0
        for segment_index in range(lines_per_kpoint):
            words = lines[line_index].split()
            line_index += 1
            values = words[1:] if segment_index == 0 else words
            for value in values:
                if eigenvalue_index >= total_eigenvalues:
                    break
                energies[ikpoint, eigenvalue_index] = float(value)
                eigenvalue_index += 1

    occupied_cutoff = fermi_level + (8.617e-5 * 300.0 * math.log(99.0))
    below_fermi = energies[energies <= occupied_cutoff]
    above_fermi = energies[energies > occupied_cutoff]
    vbm = float(np.max(below_fermi)) if below_fermi.size else fermi_level
    cbm = float(np.min(above_fermi)) if above_fermi.size else fermi_level
    return {
        "fermi_level": fermi_level,
        "vbm": vbm,
        "cbm": cbm,
        "bandgap": max(0.0, cbm - vbm),
        "nspin": nspin,
    }


def _read_pdos_energy_reference(label, work_dir, eig_path):
    bands_path = _find_optional_matching_file(label, ".bands", work_dir)
    if bands_path is not None:
        band_data = read_band_structure(bands_path)
        return {
            "fermi_level": band_data.fermi_level,
            "vbm": band_data.vbm,
            "bandgap": band_data.bandgap,
            "nspin": band_data.nspin,
            "source": bands_path,
        }

    eig_data = _read_eig_levels(eig_path)
    eig_data["source"] = eig_path
    return eig_data


def _resolve_siesta_utility(command_name):
    from NanoCore.env import siesta_util_location

    try:
        from NanoCore.env import siesta_util_pdos
    except ImportError:
        siesta_util_pdos = command_name

    utility = Path(str(siesta_util_pdos)).expanduser()
    if utility.is_absolute() or len(utility.parts) > 1:
        return utility

    utility_dir = Path(str(siesta_util_location)).expanduser()
    return utility_dir / str(siesta_util_pdos)


def _selection_from_orbital_index(orbital_index):
    tokens = str(orbital_index).split("_")
    if not tokens or any(token == "" for token in tokens) or len(tokens) > 4:
        raise ValueError(
            "Orbital selection must use atom_or_species[_n[_l[_m]]], e.g. Ba_0, C_2_1_0, or 1_2_1."
        )

    error_message = (
        "Orbital selection must use atom_or_species[_n[_l[_m]]], "
        "where n, l, and m are integers; examples: Ba_0, C_2_1_0, or 1_2_1."
    )
    selection = {
        "output": orbital_index,
        "target": tokens[0],
    }
    try:
        if len(tokens) >= 2:
            selection["n"] = int(tokens[1])
        if len(tokens) >= 3:
            selection["l"] = int(tokens[2])
        if len(tokens) >= 4:
            selection["m"] = int(tokens[3])
    except ValueError:
        raise ValueError(error_message) from None
    return selection


def _sanitize_pdos_output_token(value):
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return token.strip("._") or "selection"


def _default_pdos_output_name(selection):
    target = _sanitize_pdos_output_token(selection["target"])
    n_value = int(selection.get("n", 0))
    name_parts = ["PDOS", target, f"n{n_value}"]
    if n_value != 0:
        l_value = int(selection.get("l", -1))
        name_parts.append(f"l{l_value}")
        if l_value != -1:
            name_parts.append(f"m{int(selection.get('m', 9))}")
    return "_".join(name_parts)


def _normalize_fmpdos_selection(selection):
    if isinstance(selection, str):
        return _selection_from_orbital_index(selection)

    try:
        target = str(selection["target"]).strip()
    except (KeyError, TypeError):
        raise ValueError("PDOS selection must include a target value.") from None

    if not target:
        raise ValueError("PDOS selection target value cannot be empty.")

    normalized = {
        "target": target,
    }
    if "n" in selection and selection["n"] is not None:
        normalized["n"] = int(selection["n"])
    if "l" in selection and selection["l"] is not None:
        normalized["l"] = int(selection["l"])
    if "m" in selection and selection["m"] is not None:
        normalized["m"] = int(selection["m"])
    output_name = str(selection.get("output") or _default_pdos_output_name(normalized)).strip()
    if not output_name:
        raise ValueError("PDOS selection output value cannot be empty.")
    normalized["output"] = output_name
    return normalized


def _run_fmpdos(pdos_file, selection, executable):
    input_lines = [
        pdos_file.name,
        selection["output"],
        selection["target"],
    ]

    n_value = selection.get("n", 0)
    input_lines.append(str(n_value))
    if n_value != 0:
        l_value = selection.get("l", -1)
        input_lines.append(str(l_value))
        if l_value != -1:
            input_lines.append(str(selection.get("m", 9)))

    subprocess.run(
        [str(executable)],
        input="\n".join(input_lines) + "\n",
        text=True,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def _read_pdos_columns(filename, energy_reference, emin, emax, nspin):
    energy = []
    pdos_columns = []

    with open(filename) as pdos_file:
        for line in pdos_file:
            words = line.split()
            if not words or words[0].startswith("#"):
                continue

            shifted_energy = float(words[0]) - energy_reference
            if not (float(emin) < shifted_energy < float(emax)):
                continue

            energy.append([shifted_energy])
            row = [float(words[1])]
            if nspin > 1 and len(words) >= 3:
                row.append(-float(words[2]))
            pdos_columns.append(row)

    if not energy:
        raise ValueError(f"No PDOS data in {filename} within the requested energy window.")

    return (
        np.array(energy, dtype=float),
        np.array(pdos_columns, dtype=float),
    )


def _friendly_pdos_label(label):
    text = str(label).strip()
    if not text:
        return text

    spin_suffix = ""
    spin_match = re.search(r"\s+(spin\s+\d+)$", text, flags=re.IGNORECASE)
    if spin_match:
        spin_suffix = f" {spin_match.group(1)}"
        text = text[: spin_match.start()].strip()

    if text.lower() == "total":
        return f"Total{spin_suffix}"

    if text.startswith("PDOS_"):
        text = text[5:]

    tokens = text.split("_")
    if len(tokens) >= 4:
        try:
            target = "_".join(tokens[:-3])
            n_value = int(tokens[-3])
            l_value = int(tokens[-2])
            m_value = int(tokens[-1])
        except ValueError:
            pass
        else:
            orbital = _orbital_letter(l_value)
            return f"{target} {n_value}{orbital} m={m_value}{spin_suffix}"

    if len(tokens) < 3:
        return f"{text}{spin_suffix}"

    try:
        target = "_".join(tokens[:-2])
        n_value = int(tokens[-2])
        l_value = int(tokens[-1])
    except ValueError:
        return f"{text}{spin_suffix}"

    return f"{target} {n_value}{_orbital_letter(l_value)}{spin_suffix}"


def _orbital_letter(l_value):
    letters = {
        0: "s",
        1: "p",
        2: "d",
        3: "f",
        4: "g",
        5: "h",
    }
    return letters.get(int(l_value), f"l={int(l_value)}")


def _plot_pdos(data, labels, output_path, emin, emax):
    fig, ax = plt.subplots(figsize=(5.0, 3.0))

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    ax.tick_params(axis="x", direction="in", length=6)
    ax.tick_params(axis="y", direction="in", length=6)

    energy = data[:, 0]
    colors = ["k", "r", "tab:blue", "tab:green", "tab:orange", "tab:purple"]
    friendly_labels = [_friendly_pdos_label(label) for label in labels]
    for column_index in range(1, data.shape[1]):
        label = friendly_labels[column_index - 1] if column_index - 1 < len(friendly_labels) else None
        ax.plot(
            energy,
            data[:, column_index],
            linewidth=2,
            color=colors[(column_index - 1) % len(colors)],
            label=label,
        )

    ax.axvline(x=0.0, color="k", linestyle="--", linewidth=2.0)
    ax.set_xlim(float(emin), float(emax))
    ax.set_xticks(np.linspace(float(emin), float(emax), 5))
    ax.set_xlabel(r"$E-E_V$ (eV)", fontsize=16)
    ax.set_ylabel("DOS (states/eV)", fontsize=16)
    ax.tick_params(axis="both", labelsize=16)
    if labels:
        ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, transparent=True)
    plt.close(fig)


def generate_pdos_csv(
    pdos_path=None,
    orbital_indices=None,
    emin=-4.0,
    emax=12.0,
    output_path="PDOS.csv",
    figure_path="pdos.png",
):
    path = _find_pdos_file(pdos_path).resolve()
    work_dir = path.parent
    label = path.name[:-5] if path.name.endswith(".PDOS") else path.stem
    eig_path = _find_matching_file(label, ".EIG", work_dir)
    dos_path = _find_matching_file(label, ".DOS", work_dir)
    selected_orbitals = [_normalize_fmpdos_selection(selection) for selection in (orbital_indices or [])]

    previous_dir = Path.cwd()
    try:
        os.chdir(work_dir)
        reference = _read_pdos_energy_reference(label, Path.cwd(), Path(eig_path.name))
        fermi_level = reference["fermi_level"]
        vbm = reference["vbm"]
        nspin = reference["nspin"]
        fmpdos_executable = _resolve_siesta_utility("fmpdos")

        generated_files = []
        generated_labels = []
        for selection in selected_orbitals:
            output_file = Path(selection["output"])
            if output_file.exists():
                output_file.unlink()
            _run_fmpdos(Path(path.name), selection, fmpdos_executable)
            if not output_file.is_file():
                raise FileNotFoundError(f"fmpdos did not generate expected output file: {output_file}")
            generated_files.append(output_file)
            generated_labels.append(selection["output"])

        energy, dos_columns = _read_pdos_columns(Path(dos_path.name), vbm, emin, emax, nspin)
        data = np.hstack((energy, dos_columns))
        plot_labels = ["total"] if dos_columns.shape[1] == 1 else ["total spin 1", "total spin 2"]

        for generated_file, generated_label in zip(generated_files, generated_labels):
            _, projected_columns = _read_pdos_columns(generated_file, vbm, emin, emax, nspin)
            data = np.hstack((data, projected_columns))
            if projected_columns.shape[1] == 1:
                plot_labels.append(generated_label)
            else:
                plot_labels.extend([f"{generated_label} spin 1", f"{generated_label} spin 2"])

        np.savetxt(output_path, data, delimiter=",")
        _plot_pdos(data, plot_labels, figure_path, emin, emax)

        return {
            "csv": work_dir / output_path,
            "figure": work_dir / figure_path,
            "fermi_level": fermi_level,
            "vbm": vbm,
            "nspin": nspin,
            "pdos": path,
            "dos": dos_path,
            "reference": work_dir / reference["source"],
            "orbitals": selected_orbitals,
        }
    finally:
        os.chdir(previous_dir)


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
