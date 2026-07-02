import glob

import numpy as np
from NanoCore import s2

from .operations import siesta_eos
from .post_process import plot_band_structure, plot_pldos


BANNER = """            \\\///
           / _  _ \\       Hey, you must know what you are doing.
         (| (.)(.) |)     Otherwise you might get wrong results!
 +-----.OOOo--()--oOOO.------------------------------------------+
 |                   Python Program for SIESTA                   |
 |             py4siesta Version: 1.00 (15 July. 2020)           |
 |            Developed by RGLee    (ronggyulee@kaist.ac.kr)     |
 +-----.oooO-------------------------------------------------- --+
        (   )   Oooo.
         \\ (    (   )
          \\_)    ) /
                (_/
"""

MENU = """ ======================= K-point Sampling ========================
 1) Bulk    2) Slab    3) Analysis
 =================== Structure Optimization ======================
 4) Bulk    5) Slab    6) Sliding    7) Distance    8) Analysis
 ======================= Job Submission ==========================
 9) K-point Sampling   10) Structure Optimization
 ========================= Post-Process ==========================
 11) Band Structure    12) PLDOS
 ============================ Utility ============================
 01) Generate Geometries

 0) Quit
"""

GENERATE_GEOMETRIES_MENU = """ ======================= Generate Geometries ======================
 1) Move structure
 2) Interpolate structure

 0) Quit
"""


def _show_section(title: str) -> None:
    print(f"\n[{title}]")


def _prompt_float(message: str) -> float:
    while True:
        raw = input(message).strip()
        try:
            return float(raw)
        except ValueError:
            print("Please enter a valid number.")


def _prompt_optional_float(message: str, default: float) -> float:
    while True:
        raw = input(f"{message} [default: {default}]: ").strip()
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError:
            print("Please enter a valid number.")


def _prompt_int(message: str) -> int:
    while True:
        raw = input(message).strip()
        try:
            return int(raw)
        except ValueError:
            print("Please enter a valid integer.")


def _prompt_str(message: str) -> str:
    while True:
        value = input(message).strip()
        if value:
            return value
        print("Input cannot be empty.")


def _prompt_menu_token(message: str) -> str:
    while True:
        value = input(message).strip()
        if value:
            return value
        print("Input cannot be empty.")


def _select_existing_or_prompt(pattern: str, prompt: str):
    matches = sorted(glob.glob(pattern))
    if matches:
        file_type = pattern.replace("*", "")
        print(f"Found {file_type} files existing: {', '.join(matches)}")
        return matches[0]
    return _prompt_str(prompt)


def _prompt_choice(message: str, valid_choices) -> str:
    valid = {str(choice).strip().lower() for choice in valid_choices}
    while True:
        choice = input(message).strip().lower()
        if choice in valid:
            return choice
        print(f"Please choose one of: {', '.join(sorted(valid))}")


def _prompt_repeated_ints(message: str):
    values = []
    print("Press Enter on an empty line to finish.")
    while True:
        raw = input(message.format(index=len(values) + 1)).strip()
        if not raw:
            if values:
                return values
            print("Please enter at least one value.")
            continue
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a valid integer.")
            continue
        if value <= 0:
            print("Please enter a positive integer.")
            continue
        values.append(value)


def _prompt_direction_mask(message: str, default):
    default_mask = np.asarray(default, dtype=int).tolist()
    prompt = f"{message} [default: {' '.join(map(str, default_mask))}]: "
    while True:
        raw = input(prompt).strip()
        if not raw:
            return default_mask

        tokens = raw.replace(',', ' ').split()
        if len(tokens) != 3:
            print("Please enter exactly three values, e.g. '1 1 1' or '0 1 0'.")
            continue

        try:
            mask = [int(token) for token in tokens]
        except ValueError:
            print("Direction components must be integers 0 or 1.")
            continue

        if any(value not in (0, 1) for value in mask):
            print("Direction components must be either 0 or 1.")
            continue

        return mask


def _parse_sliding_vector(line: str) -> np.ndarray:
    normalized = line.replace(',', ' ').split()
    if len(normalized) != 2:
        raise ValueError("Each sliding vector must contain exactly two numbers.")
    return np.array([float(normalized[0]), float(normalized[1])], dtype=float)


def _prompt_sliding_vectors(mode_label: str):
    print(
        f"Enter {mode_label} sliding vectors one per line.\n"
        f"Example: {mode_label}\n"
        "Press Enter on an empty line to finish."
    )
    vectors = []
    while True:
        line = input(f"Input vector #{len(vectors) + 1}: ").strip()
        if not line:
            if vectors:
                return vectors
            print("Please enter at least one sliding vector.")
            continue
        try:
            vectors.append(_parse_sliding_vector(line))
        except ValueError as exc:
            print(exc)


def _prompt_atom_selection(struct_length: int, label: str = "atom index range to move") -> str:
    return _prompt_str(
        f"Input {label} [1-{struct_length}, e.g. 20-30]: "
    )


def _sliding_case_label(displacement_mode: str, components: np.ndarray) -> str:
    first, second = np.asarray(components, dtype=float)
    if displacement_mode == "fractional":
        return f"fa_{first:+0.4f}-fb_{second:+0.4f}"
    if displacement_mode == "absolute":
        return f"x_{first:+0.4f}-y_{second:+0.4f}"
    raise ValueError(f"Unsupported sliding mode: {displacement_mode}")


def _sliding_displacement(struct, displacement_mode: str, components: np.ndarray) -> np.ndarray:
    vector = np.asarray(components, dtype=float)
    if vector.shape != (2,):
        raise ValueError("Each sliding vector must contain exactly two in-plane components.")

    if displacement_mode == "fractional":
        cell = np.array(struct.get_cell(), dtype=float, copy=True)
        displacement = vector[0] * cell[0] + vector[1] * cell[1]
        displacement[2] = 0.0
        return displacement

    if displacement_mode == "absolute":
        return np.array([vector[0], vector[1], 0.0], dtype=float)

    raise ValueError(f"Unsupported sliding mode: {displacement_mode}")


def _prepare_sliding_cases(struct, displacement_mode: str, vectors):
    return [
        (_sliding_case_label(displacement_mode, vector), _sliding_displacement(struct, displacement_mode, vector))
        for vector in vectors
    ]


def _print_origin_structure(struct) -> None:
    print("Origin STRUCT.fdf information:")
    struct.__repr__()


def _run_generate_geometries_menu(vasp) -> None:
    print(GENERATE_GEOMETRIES_MENU)
    mode = _prompt_int("Select Generate Geometries menu number: ")

    if mode == 1:
        _show_section("Move structure")
        struct = vasp.struct
        dx = _prompt_float("Input displacement dx: ")
        dy = _prompt_float("Input displacement dy: ")
        dz = _prompt_float("Input displacement dz: ")

        struct2 = vasp.move(struct, displacement=np.array([dx, dy, dz]))
        s2.Siesta(struct2).write_struct()

    elif mode == 2:
        _show_section("Interpolate structure")
        initial_path = _prompt_str("Input initial structure path (STRUCT.fdf): ")
        final_path = _prompt_str("Input final structure path (STRUCT.fdf): ")
        division_npt = _prompt_int("Input division npt (>=2): ")
        extrapolate_npt = _prompt_int("Input extrapolate npt (0 for none): ")
        vasp.interpolate(
            initial_path=initial_path,
            final_path=final_path,
            division_npt=division_npt,
            extrapolate_npt=extrapolate_npt,
        )

    elif mode == 0:
        print("Exit Generate Geometries.")

    else:
        print("Unknown Generate Geometries menu number. Please run the program again and choose a valid option.")


def main():
    print(BANNER)
    print(MENU)
    mode_token = _prompt_menu_token("Select menu number: ")
    if mode_token == "01":
        mode = "01"
    else:
        try:
            mode = int(mode_token)
        except ValueError:
            mode = None

    vasp = None
    if mode not in {0, 11, 12, None}:
        vasp = siesta_eos()

    if mode in {4, 5, 6, 7, 8} and vasp is not None:
        _print_origin_structure(vasp.struct)
        print("\n")

    if mode == 1:
        _show_section("Bulk K-point sampling")
        kpt = _prompt_repeated_ints("Input k-point value #{index}: ")
        vasp.kpoint_sampling(kpoints=kpt)

    elif mode == 2:
        _show_section("Slab K-point sampling")
        k_values = _prompt_repeated_ints("Input kx (= ky) value #{index}: ")
        kpt = [[k, k, 1] for k in k_values]
        vasp.kpoint_sampling(sym=0, kpoints=kpt)

    elif mode == 3:
        _show_section("K-point convergence analysis")
        tolerance = _prompt_optional_float("Input total energy tolerance in eV", default=0.01)
        vasp.kpoint_analysis(tolerance=tolerance)

    elif mode == 4:
        _show_section("Bulk structure optimization")
        scale_mask = _prompt_direction_mask(
            "Input bulk expansion direction (x y z)",
            default=[1, 1, 1],
        )
        vasp.eos_bulk(scale_mask=scale_mask)

    elif mode == 5:
        _show_section("Slab structure optimization")
        include_z = _prompt_choice(
            "Include z-direction scaling? (y/n): ",
            {"y", "n"},
        )
        scale_mask = [1, 1, 1] if include_z == "y" else [1, 1, 0]
        vasp.eos_slab(scale_mask=scale_mask)

    elif mode == 6:
        _show_section("Sliding structure optimization")
        selection = _prompt_atom_selection(len(vasp.struct))
        sliding_mode = _prompt_choice(
            "Input sliding displacement mode (1: fractional, 2: absolute): ",
            {"1", "2"},
        )
        displacement_mode = "fractional" if sliding_mode == "1" else "absolute"
        mode_label = "0.25 0.50" if displacement_mode == "fractional" else "1.50 0.00"
        vectors = _prompt_sliding_vectors(mode_label)
        sliding_cases = _prepare_sliding_cases(vasp.struct, displacement_mode, vectors)
        vasp.eos_sliding(
            selection=selection,
            sliding_cases=sliding_cases,
        )

    elif mode == 7:
        _show_section("Distance structure optimization")
        selection = _prompt_atom_selection(len(vasp.struct))
        min_distance = vasp.get_distance_min(selection)
        print(f"Current minimum z distance: {min_distance:.6f}\n")

        distance_start = _prompt_float("Input distance start: ")
        distance_end = _prompt_float("Input distance end: ")
        distance_npt = _prompt_int("Input number of distance points: ")
        vasp.eos_distance(
            selection=selection,
            distance_range=np.linspace(distance_start, distance_end, distance_npt),
        )

    elif mode == 8:
        _show_section("Structure fitting")
        print("1) Bulk")
        print("2) Slab")
        print("3) Distance\n")
        select = _prompt_int("Select fitting mode: ")
        if select == 1:
            vasp.find_optimized_lattice(mode="Murnaghan")
        elif select == 2:
            vasp.find_optimized_lattice(mode="Polynomial")
        elif select == 3:
            selection = _prompt_atom_selection(len(vasp.struct), label="moving atom index range")
            vasp.find_optimized_lattice(mode="Distance", selection=selection)

    elif mode == 9:
        vasp.qsub("kpt")

    elif mode == 10:
        vasp.qsub("opt")

    elif mode == 11:
        _show_section("Plot Band Structure")
        bands_path = _select_existing_or_prompt("*.bands", "Input .bands file path: ")
        emin = _prompt_optional_float("Input minimum energy for band structure (eV)", default=-2.0)
        emax = _prompt_optional_float("Input maximum energy for band structure (eV)", default=4.0)
        result = plot_band_structure(bands_path=bands_path, emin=emin, emax=emax)
        print(f"VBM: {result['vbm']:.6f} eV")
        print(f"Band gap: {result['bandgap']:.6f} eV")
        print(f"Generated: {result['figure']}, {result['special_k']}, {result['kpath']}, {result['bands']}")

    elif mode == 12:
        _show_section("Plot PLDOS")
        pdos_path = _select_existing_or_prompt("*.PDOS", "Input .PDOS file path: ")
        emin = _prompt_optional_float("Input minimum energy for PLDOS plot (eV)", default=-4.0)
        emax = _prompt_optional_float("Input maximum energy for PLDOS plot (eV)", default=2.0)
        result = plot_pldos(pdos_path=pdos_path, emin=emin, emax=emax)
        print(f"Fermi level: {result['fermi_level']:.6f} eV")
        print(f"Generated: {result['figure']}, {result['z']}, {result['energy']}, {result['pldos']}")

    elif mode == "01":
        _run_generate_geometries_menu(vasp)

    elif mode == 0:
        print("Exit py4siesta.")

    else:
        print("Unknown menu number. Please run the program again and choose a valid option.")


if __name__ == "__main__":
    main()
