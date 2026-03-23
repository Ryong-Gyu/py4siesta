## Py4siesta Utilities

### Based on NanoCore libraries (recent version: [**OORINano**](https://github.com/yhkimlab/OORINano) from YHKLab) 


This repository provides helper scripts for SIESTA calculations, wrapping and extending the NanoCore Python tools.

### 1. Set up the environment

Prepare the NanoCore library by moving the Python library directory to your working environment.

### 2. Prepare reference calculations

Create the directory `origin` containing:

- **`input/`** – SIESTA input (`.fdf`) files  
- **`slm_siesta_run`** – Slurm job script for running SIESTA

### 3. Run Py4siesta

Execute the main script:

```bash
python py4siesta.py
```

### 4. Functionality

**1) Kpoint Sampling - Bulk**: Generate bulk calculation cases with isotropic `k x k x k` Monkhorst-Pack meshes.
**2) Kpoint Sampling - Slab**: Generate slab calculation cases with in-plane `k x k x 1` meshes for surface models.
**3) Structure Optimization - Bulk**: Create volume-scaled bulk structures for equation-of-state calculations.
**4) Structure Optimization - Slab**: Create slab structures with only the in-plane lattice vectors scaled.
**5) Structure Optimization - Sliding**: Build cases by translating selected atoms laterally to scan stacking or interface registry.
**6) Structure Optimization - Distance**: Build cases by changing the vertical distance between selected atoms and the remaining structure.
**7) Fitting**: Analyze completed EOS or distance calculations and generate an optimized structure from the fitted energy curve.
**8) Job Submission - Kpoint Sampling**: Submit all generated k-point sampling jobs with the prepared Slurm scripts.
**9) Job Submission - Structure Optimization**: Submit all generated structure optimization jobs with the prepared Slurm scripts.
**10) move**: Translate the entire structure by a user-defined displacement vector and write a new `STRUCT.fdf`.
