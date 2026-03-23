import numpy as np
from NanoCore import s2

from operations import siesta_eos


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

MENU = """ ======================= Kpoint Sampling =========================
 1) Bulk    2) Slab
 =================== Structure Optimization ======================
 3) Bulk    4) Slab    5) Sliding    6) Distance    7) Fitting
 ======================= Job Submission ==========================
 8) Kpoint Sampling    9) Structure Optimization
 ============================ Utility ============================
10) move

 0) Quit
"""


def _prompt_float(message: str) -> float:
    return float(input(message))


def _prompt_int(message: str) -> int:
    return int(input(message))


def _prompt_str(message: str) -> str:
    return input(message).strip()


def _prompt_choice(message: str, valid_choices) -> str:
    valid = {str(choice).strip().lower() for choice in valid_choices}
    while True:
        choice = input(message).strip().lower()
        if choice in valid:
            return choice
        print(f"Please choose one of: {', '.join(sorted(valid))}")


def _parse_sliding_vector(line: str) -> np.ndarray:
    normalized = line.replace(',', ' ').split()
    if len(normalized) != 2:
        raise ValueError("Each sliding vector must contain exactly two numbers.")
    return np.array([float(normalized[0]), float(normalized[1])], dtype=float)


def _prompt_sliding_vectors(mode_label: str):
    print(
        f"Enter {mode_label} sliding vectors one per line in the form 'v1 v2'.\n"
        "Press Enter on an empty line to finish."
    )
    vectors = []
    while True:
        line = input(f"{mode_label} vector #{len(vectors) + 1}: ").strip()
        if not line:
            if vectors:
                return vectors
            print("Please enter at least one sliding vector.")
            continue
        try:
            vectors.append(_parse_sliding_vector(line))
        except ValueError as exc:
            print(exc)


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


def main():
    print(BANNER)
    print(MENU)
    mode = _prompt_int("")

    vasp = siesta_eos()

    if mode in {3, 4, 5, 6, 7}:
        _print_origin_structure(vasp.struct)
        print("\n")

    if mode == 1:
        kpt = []
        while True:
            k = _prompt_int("Type number of k points for calculation (0: quit): ")
            if k == 0:
                break
            elif k > 0:
                kpt.append(k)
        vasp.kpoint_sampling(kpoints=kpt)

    elif mode == 2:
        kpt = []
        while True:
            k = _prompt_int("Type number of kx (=ky) for calculation (0: quit): ")
            if k == 0:
                break
            elif k > 0:
                kpt.append([k, k, 1])
        vasp.kpoint_sampling(sym=0, kpoints=kpt)

    elif mode == 3:
        vasp.eos_bulk()

    elif mode == 4:
        vasp.eos_slab()

    elif mode == 5:
        selection = _prompt_str("Moving atom index range (e.g. 21-30 (from 21 to 30 atoms)): ")
        sliding_mode = _prompt_choice(
            "Sliding displacement mode (1: fractional in a/b, 2: absolute in Ang): ",
            {"1", "2"},
        )
        displacement_mode = "fractional" if sliding_mode == "1" else "absolute"
        mode_label = "f_a f_b" if displacement_mode == "fractional" else "x y"
        vectors = _prompt_sliding_vectors(mode_label)
        sliding_cases = _prepare_sliding_cases(vasp.struct, displacement_mode, vectors)
        vasp.eos_sliding(
            selection=selection,
            sliding_cases=sliding_cases,
        )

    elif mode == 6:
        min_distance = vasp.get_distance_min(selection)
        print(f"Current minimum z distance (in Ang): {min_distance:.6f}\n")

        selection = _prompt_str("Moving atom index range (e.g. 21-30 (from 21 to 30 atoms)): ")
        distance_start = _prompt_float("Distance start: ")
        distance_end = _prompt_float("Distance end: ")
        distance_npt = _prompt_int("Distance npt: ")
        vasp.eos_distance(
            selection=selection,
            distance_range=np.linspace(distance_start, distance_end, distance_npt),
        )

    elif mode == 7:
        print("1) Bulk\n")
        print("2) Slab\n")
        print("3) Distance\n")
        select = _prompt_int(": ")
        if select == 1:
            vasp.find_optimized_lattice(mode="Murnaghan")
        elif select == 2:
            vasp.find_optimized_lattice(mode="Polynomial")
        elif select == 3:
            selection = _prompt_str("Moving atom index range (e.g. 20-30): ")
            vasp.find_optimized_lattice(mode="Distance", selection=selection)

    elif mode == 8:
        vasp.qsub("kpt")

    elif mode == 9:
        vasp.qsub("opt")

    elif mode == 10:
        struct = vasp.struct
        while True:
            dx = _prompt_float("displacement dx: ")
            dy = _prompt_float("displacement dy: ")
            dz = _prompt_float("displacement dz: ")
            break

        struct2 = vasp.move(struct, displacement=np.array([dx, dy, dz]))
        s2.Siesta(struct2).write_struct()

    elif mode == 0:
        pass


if __name__ == "__main__":
    main()
