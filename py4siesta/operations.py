import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import fminbound, leastsq

from NanoCore import *

from .utils import copy_contents, last_matching_line, working_dir


class SiestaContext:
    def __init__(self):
        self.root = Path.cwd()
        self.origin_dir = self.root / "origin"
        self.struct = s2.read_fdf(self.origin_dir / "input" / "STRUCT.fdf")


class BaseOperation:
    base_dirname = ""
    output_name = "STRUCT.fdf"

    def __init__(self, context: SiestaContext):
        self.context = context

    @property
    def base_dir(self):
        return self.context.root / self.base_dirname

    @staticmethod
    def _normalize_scale_mask(scale_mask, default):
        if scale_mask is None:
            scale_mask = default

        values = np.asarray(scale_mask, dtype=int).reshape(-1)
        if values.size != 3:
            raise ValueError("Scale direction must contain exactly three components.")
        if not np.all(np.isin(values, [0, 1])):
            raise ValueError("Scale direction components must be either 0 or 1.")
        if not np.any(values):
            raise ValueError("At least one scale direction must be enabled.")
        return values.astype(float)

    @classmethod
    def _scale_case(cls, ratio, scale_mask, default_mask):
        mask = cls._normalize_scale_mask(scale_mask, default_mask)
        return {"ratio": float(ratio), "mask": mask}

    @staticmethod
    def _scaling_factors(ratio, mask):
        mask_array = np.asarray(mask, dtype=float)
        return 1.0 + (float(ratio) - 1.0) * mask_array

    @classmethod
    def _scaled_structure(cls, struct, ratio, scale_mask, default_mask):
        mask = cls._normalize_scale_mask(scale_mask, default_mask)
        factors = cls._scaling_factors(ratio, mask)
        scaled_struct = struct.copy()

        cell = np.array(struct.get_cell(), dtype=float, copy=True)
        scaled_struct.set_cell(cell * factors[np.newaxis, :])

        for atom, scaled_atom in zip(struct, scaled_struct._atoms):
            position = np.array(atom.get_position(), dtype=float, copy=True)
            scaled_atom.set_position(Vector(position * factors))

        return scaled_struct

    def write_metadata(self, **kwargs):
        return None

    @staticmethod
    def _print_structure_summary(struct):
        print("Origin cell:")
        print(np.array(struct.get_cell(), copy=True))
        print("Origin atomic geometry:")
        for atom in struct:
            position = np.array(atom.get_position(), dtype=float)
            print(
                f"{atom.get_serial():4d} {atom.get_symbol():>2s} "
                f"{position[0]:12.6f} {position[1]:12.6f} {position[2]:12.6f}"
            )

    @staticmethod
    def _selected_serials(struct, selection):
        selected = struct.copy()
        selected.select_atmnbs(selection)
        return list(selected._selected)

    @classmethod
    def _validate_selection(cls, struct, selection):
        selected_serials = cls._selected_serials(struct, selection)
        all_serials = set(struct.get_serials())
        remaining_serials = sorted(all_serials - set(selected_serials))
        if not selected_serials:
            raise ValueError("No atoms were selected.")
        if not remaining_serials:
            raise ValueError("Selection must leave at least one fixed atom.")
        return selected_serials, remaining_serials

    @staticmethod
    def _translate_selected(struct, selected_serials, displacement):
        translated = struct.copy()
        translated._selected = list(selected_serials)
        translated.translate(*np.asarray(displacement, dtype=float))
        return translated

    @staticmethod
    def _minimum_z_distance(struct, selected_serials, remaining_serials):
        selected_z = [float(struct._atoms[serial - 1].get_position()[2]) for serial in selected_serials]
        remaining_z = [float(struct._atoms[serial - 1].get_position()[2]) for serial in remaining_serials]
        return min(abs(z_sel - z_fix) for z_sel in selected_z for z_fix in remaining_z)

    def prepare_base_dir(self):
        base_dir = self.base_dir
        if base_dir.exists():
            shutil.rmtree(base_dir)
        base_dir.mkdir()
        return base_dir

    def iter_cases(self, case_parameters):
        for index, case_parameter in enumerate(case_parameters, start=1):
            yield index, case_parameter

    def prepare_case_dir(self, case_name):
        return working_dir(Path(case_name))

    def copy_origin(self):
        copy_contents(self.context.origin_dir, Path.cwd())

    def finalize_case(self, case_input):
        destination = Path("input")
        if self.output_name == "STRUCT.fdf":
            s2.Siesta(case_input).write_struct()
        elif self.output_name == "KPT.fdf":
            simulation = s2.Siesta(self.context.struct)
            simulation.set_option("kgrid", list(case_input))
            simulation.write_kpt()
        else:
            raise ValueError(f"Unsupported output name: {self.output_name}")

        shutil.move(self.output_name, destination / self.output_name)

    def case_parameters(self, **kwargs):
        raise NotImplementedError

    def case_name(self, index, case_parameter):
        raise NotImplementedError

    def build_case_input(self, case_parameter):
        raise NotImplementedError

    def run(self, **kwargs):
        self.prepare_base_dir()

        with working_dir(self.base_dir):
            self.write_metadata(**kwargs)
            for index, case_parameter in self.iter_cases(self.case_parameters(**kwargs)):
                with self.prepare_case_dir(self.case_name(index, case_parameter)):
                    self.copy_origin()
                    self.finalize_case(self.build_case_input(case_parameter))


class KPointSamplingOperation(BaseOperation):
    base_dirname = "01.kpoint_sampling"
    output_name = "KPT.fdf"

    def case_parameters(self, sym=1, kpoints=None):
        if kpoints is None:
            kpoints = [1, 1, 1]

        if sym == 1:
            return [[k, k, k] for k in kpoints]
        return [list(k) for k in kpoints]

    def case_name(self, index, case_parameter):
        return f"{case_parameter[0]}+{case_parameter[1]}+{case_parameter[2]}"

    def build_case_input(self, case_parameter):
        return case_parameter


class BulkEosOperation(BaseOperation):
    base_dirname = "02.volume_eos"
    default_scale_mask = np.array([1.0, 1.0, 1.0])

    def write_metadata(self, ratio_range=None, scale_mask=None):
        metadata = {
            "mode": "bulk",
            "scale_mask": self._normalize_scale_mask(scale_mask, self.default_scale_mask).astype(int).tolist(),
        }
        Path("eos_config.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def case_parameters(self, ratio_range=None, scale_mask=None):
        if ratio_range is None:
            ratio_range = np.linspace(0.99, 1.01, 11)
        return [self._scale_case(ratio, scale_mask, self.default_scale_mask) for ratio in ratio_range]

    def case_name(self, index, case_parameter):
        return f"{index:02d}-{case_parameter['ratio']:4.3f}"

    def build_case_input(self, case_parameter):
        return self._scaled_structure(
            self.context.struct,
            ratio=case_parameter["ratio"],
            scale_mask=case_parameter["mask"],
            default_mask=self.default_scale_mask,
        )


class SlabEosOperation(BaseOperation):
    base_dirname = "02.slab_eos"
    default_scale_mask = np.array([1.0, 1.0, 0.0])

    def write_metadata(self, ratio_range=None, scale_mask=None):
        metadata = {
            "mode": "slab",
            "scale_mask": self._normalize_scale_mask(scale_mask, self.default_scale_mask).astype(int).tolist(),
        }
        Path("eos_config.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def case_parameters(self, ratio_range=None, scale_mask=None):
        if ratio_range is None:
            ratio_range = np.linspace(0.98, 1.02, 11)
        return [self._scale_case(ratio, scale_mask, self.default_scale_mask) for ratio in ratio_range]

    def case_name(self, index, case_parameter):
        return f"{index:02d}-{case_parameter['ratio']:4.3f}"

    def build_case_input(self, case_parameter):
        return self._scaled_structure(
            self.context.struct,
            ratio=case_parameter["ratio"],
            scale_mask=case_parameter["mask"],
            default_mask=self.default_scale_mask,
        )


class SlidingOperation(BaseOperation):
    base_dirname = "02.sliding"

    def case_parameters(self, selection, sliding_cases):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, _ = self._validate_selection(struct, selection)
        return [
            (selected_serials, str(case_label), np.asarray(displacement, dtype=float))
            for case_label, displacement in sliding_cases
        ]

    def case_name(self, index, case_parameter):
        _, case_label, _ = case_parameter
        return f"{index:02d}-{case_label}"

    def build_case_input(self, case_parameter):
        selected_serials, _, displacement = case_parameter
        return self._translate_selected(self.context.struct, selected_serials, displacement)


class DistanceOperation(BaseOperation):
    base_dirname = "02.distance"

    def current_distance(self, selection):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, remaining_serials = self._validate_selection(struct, selection)
        return self._minimum_z_distance(struct, selected_serials, remaining_serials)

    def case_parameters(self, selection, distance_range):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, remaining_serials = self._validate_selection(struct, selection)
        current_distance = self._minimum_z_distance(struct, selected_serials, remaining_serials)
        return [
            (selected_serials, remaining_serials, current_distance, float(target_distance))
            for target_distance in distance_range
        ]

    def case_name(self, index, case_parameter):
        _, _, _, target_distance = case_parameter
        return f"{index:02d}-distance_{target_distance:0.4f}"

    def build_case_input(self, case_parameter):
        selected_serials, _, current_distance, target_distance = case_parameter
        dz = target_distance - current_distance
        return self._translate_selected(self.context.struct, selected_serials, [0.0, 0.0, dz])


class FitOptimizedStructureOperation:
    def __init__(self, context: SiestaContext):
        self.context = context

    @staticmethod
    def _distance_from_case_name(case_name):
        return float(case_name.split("distance_")[-1])

    @staticmethod
    def _load_scale_mask(base_dir, default_mask):
        config_path = base_dir / "eos_config.json"
        if not config_path.is_file():
            return np.asarray(default_mask, dtype=float)

        metadata = json.loads(config_path.read_text())
        return BaseOperation._normalize_scale_mask(metadata.get("scale_mask"), default_mask)

    def run(self, mode="Murnaghan", selection=None):
        struct = self.context.struct.copy()
        atoms = struct._atoms
        cell = np.array(struct.get_cell(), copy=True)
        init_volume = abs(np.dot(cell[2], np.cross(cell[0], cell[1])))
        init_lattice = np.sqrt(np.dot(cell[0], cell[0]))

        if mode == "Murnaghan":
            base_dir = self.context.root / "02.volume_eos"
            scale_mask = self._load_scale_mask(base_dir, BulkEosOperation.default_scale_mask)
        elif mode == "Polynomial":
            base_dir = self.context.root / "02.slab_eos"
            scale_mask = self._load_scale_mask(base_dir, SlabEosOperation.default_scale_mask)
        elif mode == "Distance":
            base_dir = self.context.root / "02.distance"
            scale_mask = None
        else:
            base_dir = self.context.root
            scale_mask = None

        if not base_dir.exists():
            raise FileNotFoundError(f"{base_dir} does not exist")

        energy = []
        volume = []
        lattice = []

        optimized_dir = base_dir / "optimized_structure"
        with working_dir(base_dir):
            if optimized_dir.exists():
                shutil.rmtree(optimized_dir)

            for path in sorted(base_dir.iterdir()):
                if not path.is_dir() or path.name == "optimized_structure":
                    continue

                out_dir = path / "OUT"
                stdout_path = out_dir / "stdout.txt"

                if not stdout_path.is_file():
                    continue

                print(path.name)
                energy_line = last_matching_line(stdout_path, "siesta:         Total =")
                volume_line = last_matching_line(stdout_path, "outcell: Cell volume")
                lattice_line = last_matching_line(stdout_path, "outcell: Cell vector modules")

                if not energy_line or not volume_line:
                    continue

                energy.append(float(energy_line.split()[-1]))
                volume.append(float(volume_line.split()[-1]))

                if mode == "Distance":
                    lattice.append(self._distance_from_case_name(path.name))
                elif lattice_line:
                    lattice.append(float(lattice_line.split()[-3]))

        energy = np.array(energy, dtype=float)
        volume = np.array(volume, dtype=float)
        lattice = np.array(lattice, dtype=float)

        a, b, c = plt.polyfit(volume, energy, 2)
        coeff_poly_4nd = plt.polyfit(lattice, energy, 4)

        v0 = -b / (2 * a)
        e0 = a * v0 ** 2 + b * v0 + c
        b0 = 2 * a * v0
        bP = 4
        x0 = [e0, b0, bP, v0]

        def murnaghan(parameters, vol):
            e0_local = parameters[0]
            b0_local = parameters[1]
            bp_local = parameters[2]
            v0_local = parameters[3]
            return e0_local + b0_local * vol / bp_local * (((v0_local / vol) ** bp_local) / (bp_local - 1) + 1) - v0_local * b0_local / (bp_local - 1.0)

        def polynomial(parameters, x):
            a_local = parameters[0]
            b_local = parameters[1]
            c_local = parameters[2]
            d_local = parameters[3]
            e_local = parameters[4]
            return a_local * x ** 4 + b_local * x ** 3 + c_local * x ** 2 + d_local * x + e_local

        func_poly_4nd = lambda x: polynomial(coeff_poly_4nd, x)

        def loss_function(parameters, y, x):
            return y - murnaghan(parameters, x)

        if mode == "Murnaghan":
            vfit = np.linspace(min(volume), max(volume), 100)
            opt_coeff, ier = leastsq(loss_function, x0, args=(energy, volume))
            opt_volume = opt_coeff[3]
            opt_func = murnaghan(opt_coeff, vfit)

            active_directions = int(np.count_nonzero(scale_mask))
            ratio = (opt_volume / init_volume) ** (1 / active_directions)
            factors = BaseOperation._scaling_factors(ratio, scale_mask)

            for iatom in range(len(atoms)):
                pos = factors * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            struct._cell = cell * factors[np.newaxis, :]
            plt.plot(volume, energy, "ro")

        elif mode == "Polynomial":
            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_lattice = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_func = polynomial(coeff_poly_4nd, vfit)

            ratio = opt_lattice / init_lattice
            factors = BaseOperation._scaling_factors(ratio, scale_mask)

            for iatom in range(len(atoms)):
                pos = factors * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            struct._cell = cell * factors[np.newaxis, :]
            plt.plot(lattice, energy, "ro")

        elif mode == "Distance":
            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_distance = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_func = polynomial(coeff_poly_4nd, vfit)
            distance_operation = DistanceOperation(self.context)
            if selection is None:
                raise ValueError("Distance fitting requires a moving atom selection.")
            selected_serials, remaining_serials = distance_operation._validate_selection(
                self.context.struct,
                selection,
            )
            current_distance = distance_operation._minimum_z_distance(
                self.context.struct,
                selected_serials,
                remaining_serials,
            )
            struct = distance_operation._translate_selected(
                self.context.struct,
                selected_serials,
                [0.0, 0.0, opt_distance - current_distance],
            )
            struct._cell = cell.copy()
            plt.plot(lattice, energy, "ro")

        with working_dir(base_dir):
            plt.plot(vfit, opt_func)
            plt.savefig("eos_fitting.png")

            copy_contents(self.context.origin_dir, optimized_dir)
            with working_dir(optimized_dir):
                s2.Siesta(struct).write_struct()
                shutil.move("STRUCT.fdf", Path("input") / "STRUCT.fdf")


class JobSubmissionOperation:
    def __init__(self, context: SiestaContext):
        self.context = context

    def run(self, mode):
        if mode == "kpt":
            targets = sorted(self.context.root.glob("01.*"))
        elif mode == "opt":
            targets = sorted(self.context.root.glob("02.*"))
        else:
            targets = []

        if not targets:
            return

        base_dir = targets[0]
        with working_dir(base_dir):
            for subdir in sorted(Path.cwd().iterdir()):
                if not subdir.is_dir():
                    continue
                with working_dir(subdir):
                    for script in sorted(Path.cwd().glob("slm_*")):
                        subprocess.run(["sbatch", str(script)], check=True)


class MoveStructureOperation:
    def __init__(self, context: SiestaContext):
        self.context = context

    def run(self, struct, displacement=np.array([0, 0, 0])):
        struct2 = struct.copy()
        struct2.select_all()
        struct2.translate(*np.asarray(displacement, dtype=float))
        return struct2


class siesta_eos:
    def __init__(self):
        self.context = SiestaContext()
        self.root = self.context.root
        self.origin_dir = self.context.origin_dir
        self.struct = self.context.struct
        self._kpoint_sampling = KPointSamplingOperation(self.context)
        self._bulk_eos = BulkEosOperation(self.context)
        self._slab_eos = SlabEosOperation(self.context)
        self._sliding = SlidingOperation(self.context)
        self._distance = DistanceOperation(self.context)
        self._fit_optimized_structure = FitOptimizedStructureOperation(self.context)
        self._job_submission = JobSubmissionOperation(self.context)
        self._move_structure = MoveStructureOperation(self.context)

    def kpoint_sampling(self, sym=1, kpoints=None):
        return self._kpoint_sampling.run(sym=sym, kpoints=kpoints)

    def eos_bulk(self, ratio_range=None, scale_mask=None):
        return self._bulk_eos.run(ratio_range=ratio_range, scale_mask=scale_mask)

    def eos_slab(self, ratio_range=None, scale_mask=None):
        return self._slab_eos.run(ratio_range=ratio_range, scale_mask=scale_mask)

    def eos_sliding(self, selection, sliding_cases):
        return self._sliding.run(selection=selection, sliding_cases=sliding_cases)

    def get_distance_min(self, selection):
        return self._distance.current_distance(selection)

    def eos_distance(self, selection, distance_range):
        return self._distance.run(selection=selection, distance_range=distance_range)

    def find_optimized_lattice(self, mode="Murnaghan", selection=None):
        return self._fit_optimized_structure.run(mode=mode, selection=selection)

    def qsub(self, mode):
        return self._job_submission.run(mode)

    def move(self, struct, displacement=np.array([0, 0, 0])):
        return self._move_structure.run(struct, displacement=displacement)
