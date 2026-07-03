# py4siesta

`py4siesta` is a numbered-menu command-line tool for preparing, organizing, submitting, and post-processing SIESTA calculation cases.

It is designed around one workflow:

1. Prepare a reference SIESTA input in an `origin/` directory.
2. Start the CLI.
3. Choose a numbered menu option.
4. Let `py4siesta` generate or analyze the calculation cases.

## Important Notice

This project includes the bundled `NanoCore` code developed by Dr. Hu Sung Kim, KISTI, and related original work from **OORINano / YHKLab**.

- Keep the original attribution and copyright context for `NanoCore` and OORINano.
- The local bundled `NanoCore/` package is used directly from this repository.
- You do not need to download `NanoCore` from an external source.

## Quick Start

Clone and enter the repository:

```bash
git clone https://github.com/YHKlab-RGLee/py4siesta.git
cd py4siesta
```

Install the package:

```bash
pip install .
```

For development, use editable install instead:

```bash
pip install -e .
```

Run the CLI:

```bash
py4siesta
```

Or run it directly from the repository root:

```bash
python -m py4siesta
```

## Requirements

Runtime Python dependencies are declared in `pyproject.toml`:

- `numpy`
- `scipy`
- `matplotlib`

The repository also installs the bundled `NanoCore` package.

## Configure SIESTA Paths

The bundled `NanoCore` module reads local executable paths from `NanoCore/env.py`. Update this file before using features that call SIESTA utilities.

Common values to edit:

```python
siesta_calculator = "/path/to/siesta"
siesta_psf_location = "/path/to/psf"
siesta_util_location = "/path/to/siesta/utils"
siesta_util_band = "gnubands"
siesta_util_dos = "Eig2DOS"
siesta_util_pdos = "fmpdos"
siesta_util_rho = "rho2xsf"
siesta_util_vh = "macroave"
siesta_util_pldos = "/path/to/PyProjection.py"
```

For PDOS extraction, make sure `siesta_util_location` points to the directory containing `fmpdos`, or set `siesta_util_pdos` to the full executable path.

## Working Directory Setup

Before running case-generation menus, create an `origin/` directory in your working directory.

Expected layout:

```text
origin/
├── input/
│   └── STRUCT.fdf
└── slm_siesta_run
```

Required inputs:

| Path | Purpose |
| --- | --- |
| `origin/input/STRUCT.fdf` | Reference SIESTA structure/input file copied into generated cases. |
| `origin/slm_*` | One or more job scripts copied into generated cases. |

## Menu Reference

Start the program, then select a menu number:

```bash
py4siesta
```

### K-point Sampling

| Menu | Feature | Purpose |
| --- | --- | --- |
| `1` | Bulk | Generate bulk cases using isotropic `k x k x k` Monkhorst-Pack meshes. |
| `2` | Slab | Generate slab cases using in-plane `k x k x 1` meshes. |
| `3` | Analysis | Run k-point convergence analysis from completed k-point calculations. |

### Structure Optimization

| Menu | Feature | Purpose |
| --- | --- | --- |
| `4` | Bulk | Create bulk equation-of-state structures and choose scaled directions (`x y z`). |
| `5` | Slab | Create slab equation-of-state structures with in-plane-focused scaling. |
| `6` | Sliding | Translate selected atoms laterally to scan stacking or interface registry. |
| `7` | Distance | Change the vertical distance between selected atoms and the remaining structure. |
| `8` | Analysis | Analyze completed EOS or distance calculations and generate an optimized structure. |

### Job Submission

| Menu | Feature | Purpose |
| --- | --- | --- |
| `9` | K-point Sampling | Submit generated k-point sampling jobs. |
| `10` | Structure Optimization | Submit generated structure-optimization jobs. |

### Post-Process

| Menu | Feature | Inputs | Outputs |
| --- | --- | --- | --- |
| `11` | Plot Band Structure | A `.bands` file and energy window. | `band.png`, plus band-path helper files. |
| `12` | Generate PDOS CSV | SIESTA `.PDOS`, matching `.EIG` and `.DOS`, energy window, and `fmpdos` atom/orbital selections. | `PDOS.csv` and `pdos.png`. |
| `13` | Plot PLDOS | A `.PDOS` file, matching `.xyz` and `.EIG`, and energy window. | `pldos.png`, plus PLDOS helper files. |

Notes:

- Menu `11` automatically selects a local `.bands` file when one is present.
- Menu `12` uses a matching `<label>.bands` file for `E_v` when available; otherwise it estimates `E_v` from `<label>.EIG`.
- Menu `12` removes existing same-named `fmpdos` output files before regenerating PDOS data.
- Menu `13` automatically selects a local `.PDOS` file when one is present.

### Utility

| Menu | Feature | Purpose |
| --- | --- | --- |
| `01` | Generate Geometries | Open the geometry-generation submenu. |

Geometry-generation submenu:

| Submenu | Feature | Purpose |
| --- | --- | --- |
| `1` | Move structure | Translate the full structure by a user-defined displacement vector and write a new `STRUCT.fdf`. |
| `2` | Interpolate structure | Generate interpolated structures between two `STRUCT.fdf` files, with optional extrapolated points. |

## Running Modes

Use one of these commands depending on how you installed or checked out the project:

```bash
# From an installed package
py4siesta

# From the repository root
python -m py4siesta

# Direct script usage
python {py4siesta_path}/run.py
```

## Agent JSON Tool

For AI agents or scripts that need non-interactive access while preserving the numbered menu CLI, use:

```bash
py4siesta-tool <command> [options]
```

The command prints one JSON object. Successful commands include `"ok": true`; failures return a nonzero exit status and a JSON error object.

Examples:

```bash
py4siesta-tool kpoint-bulk --kpoints 2 4 6
py4siesta-tool eos-sliding --selection 20-30 --mode fractional --vector 0.25 0.50
py4siesta-tool band --bands-path Solid.bands --emin -2 --emax 4
py4siesta-tool pdos --pdos-path MgO.PDOS --orbital Mg_0 --emin -4 --emax 12
py4siesta-tool pldos --pdos-path Device.PDOS --emin -4 --emax 2
```

Run `py4siesta-tool --help` or `py4siesta-tool <command> --help` for the full command list and required inputs.

## Bundled NanoCore

This project is based on `NanoCore`, an older version of the ongoing [`OORINano`](https://github.com/yhkimlab/OORINano) project. `NanoCore` is imported locally from this repository.

That means:

- `python -m py4siesta` works from this repository because `NanoCore/` is present locally.
- `pip install .` or `pip install -e .` installs both `py4siesta` and the bundled `NanoCore` package defined in `pyproject.toml`.
- There is no need to fetch `NanoCore` from an external site or package index.
