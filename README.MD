# py4siesta

`py4siesta` is a CLI utility for preparing and organizing SIESTA calculation cases.

## Important notice

This project is based on the bundled `NanoCore` code (developed by Dr. Hu Sung Kim, KISTI) included in this repository and related original work from **OORINano / YHKLab**.

- Do **not** remove the original attribution/copyright context for `NanoCore` and OORINano.
- This repository uses the local bundled `NanoCore` package; it does **not** require downloading `NanoCore` from an external source.

## Installation

### Clone the repository

```bash
git clone https://github.com/YHKlab-RGLee/py4siesta.git
cd py4siesta
```

### Install

Standard install:

```bash
pip install .
```

Editable install:

```bash
pip install -e .
```

## Dependencies

The package declares these runtime dependencies:

- `numpy`
- `scipy`
- `matplotlib`

The repository also includes the private/bundled `NanoCore` module used by the CLI.

### About `NanoCore`

This project is based on `NanoCore` which is old version of on-going project [`OORINano`](https://github.com/yhkimlab/OORINano). `NanoCore` is imported locally from this repository.

That means:

- `python -m py4siesta` works from this repository because the bundled `NanoCore/` package is present locally
- `pip install .` or `pip install -e .` from this repository installs both `py4siesta` and the bundled `NanoCore` package defined in `pyproject.toml`
- there is no need to fetch `NanoCore` from an external site or package index

## Running the program

Run from the repository root as a module:

```bash
python -m py4siesta
```

Or, after installation, run the console command:

```bash
py4siesta
```

## Required working directory layout

Before running the program, create an `origin` directory in your working directory.

Example:

```text
origin/
├── input/
│   └── STRUCT.fdf
└── slm_siesta_run
```

Required files:

- `origin/input/STRUCT.fdf`
- one or more job scripts in `origin` matching `slm_*`

The tool copies this reference input and creates calculation case directories from it.

## Features


### 1. K-point Sampling
- **Menu 1 - Bulk**
  Generate bulk calculation cases using isotropic `k x k x k` Monkhorst-Pack meshes.
- **Menu 2 - Slab**
  Generate slab calculation cases using in-plane `k x k x 1` meshes.
- **Menu 3 - Analysis**
  Perform k-point convergence test from the completed k-point sampling calculations.

### 2. Structure Optimization
- **Menu 4 - Bulk**
   Create bulk EOS structures while choosing which lattice directions (`x y z`) are scaled.
- **Menu 5 - Slab**
   Create slab EOS structures with in-plane-focused scaling.
- **Menu 6 - Sliding**
   Build cases by translating selected atoms laterally to scan stacking or interface registry.
- **Menu 7 - Distance**
   Build cases by changing the vertical distance between selected atoms and the remaining structure.
- **Menu 8 - Analysis**
   Analyze completed EOS or distance calculations and generate an optimized structure from the fitted energy curve.

### 3. Job Submission
- **Menu 9 - K-point Sampling**
   Submit generated k-point sampling jobs.
- **Menu 10 - Structure Optimization**
   Submit generated structure-optimization jobs.

### 4. Extra
- **Menu 01 - Generate Geometries**
   Open the geometry-generation submenu.
- **Submenu 1 - Move structure**
   Translate the entire structure by a user-defined displacement vector and write a new `STRUCT.fdf`.
- **Submenu 2 - Interpolate structure**
   Generate interpolated structures between two `STRUCT.fdf` files, with optional extrapolated points.

## Quick usage summary

Development usage:

```bash
python -m py4siesta
```

Installed usage:  

```bash
py4siesta
```

Direct usage:  

```bash
python {py4siesta_path}/run.py
```
