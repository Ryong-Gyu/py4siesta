# py4siesta

`py4siesta` is a CLI utility for preparing and organizing SIESTA calculation cases.

## Important notice

This project is based on the bundled `NanoCore` code included in this repository and related original work from **OORINano / YHKLab**.

- Do **not** remove the original attribution/copyright context for `NanoCore` and OORINano.
- This repository uses the local bundled `NanoCore` package; it does **not** require downloading `NanoCore` from an external source.

## Installation

### Clone the repository

```bash
git clone <repo-url>
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

`NanoCore` is imported locally from this repository.

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

1. **K-point Sampling - Bulk**  
   Generate bulk calculation cases using isotropic `k x k x k` Monkhorst-Pack meshes.
2. **K-point Sampling - Slab**  
   Generate slab calculation cases using in-plane `k x k x 1` meshes.
3. **Structure Optimization - Bulk**  
   Create bulk EOS structures while choosing which lattice directions (`x y z`) are scaled.
4. **Structure Optimization - Slab**  
   Create slab EOS structures with in-plane-focused scaling.
5. **Structure Optimization - Sliding**  
   Build cases by translating selected atoms laterally to scan stacking or interface registry.
6. **Structure Optimization - Distance**  
   Build cases by changing the vertical distance between selected atoms and the remaining structure.
7. **Fitting**  
   Analyze completed EOS or distance calculations and generate an optimized structure from the fitted energy curve.
8. **Job Submission - K-point Sampling**  
   Submit generated k-point sampling jobs.
9. **Job Submission - Structure Optimization**  
   Submit generated structure-optimization jobs.
10. **Move structure**  
    Translate the entire structure by a user-defined displacement vector and write a new `STRUCT.fdf`.

## Quick usage summary

Development usage:

```bash
python -m py4siesta
```

Installed usage:

```bash
py4siesta
```
