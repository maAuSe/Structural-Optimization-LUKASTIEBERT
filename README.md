# Structural Optimization - MTOP Assignment

This repository contains the assignment material, lecture code, references, and MATLAB workspace for the course *Structural Optimization*. The assignment is to assess the computational efficiency of multiresolution topology optimization (MTOP) for 2D minimum-compliance topology optimization of the half MBB beam.

## Assignment Snapshot

- Topic: multiresolution topology optimization (MTOP) based on Nguyen et al. (2010).
- Main comparison: classical density-based topology optimization versus MTOP.
- Structural model: 2D half MBB beam, minimum compliance, volume constraint.
- Report deadline: 15 May 2026.
- Presentation date: 22 May 2026.
- Submission package: report PDF plus a ZIP archive containing all MATLAB files and a brief `readme.txt`.

The assignment asks for four main studies:

1. Run the classical MBB beam code with optimality criteria (OC), `600 x 200` finite elements, volume fraction `0.5`, penalization `3`, and sensitivity filtering with filter radius `24` elements. Plot compliance and volume fraction against iteration number and computation time.
2. Implement MTOP with a `120 x 40` finite element mesh and a `5 x 5` density mesh inside each finite element. Compare the resulting `600 x 200` density design with the classical result.
3. Repeat the classical-versus-MTOP comparison with density filtering.
4. Switch to MMA and investigate MTOP performance for sensitivity filtering, density filtering, and Heaviside projection with different projection thresholds.

## Repository Layout

| Path | Contents |
| --- | --- |
| `assignment_summary.md` | Parsed assignment brief, deadlines, report requirements, and problem statement. |
| `course_slides/` | Course slides used as theory references. |
| `relevant_papers/` | MTOP paper and related references. |
| `structopt_exercises_lectures/` | Lecture exercise MATLAB scripts. Chapter 9 contains the topology optimization starting point. |
| `structopt/` | Optimization utilities, including `mma.m`, `gcmma.m`, and finite-difference check plotting helpers. |
| `stabil/` | STABiL finite element toolbox used by other course examples. |
| `matlab/` | Assignment workspace for new code, figures, and generated results. |
| `report/` | LaTeX report workspace with a structured `main.tex` draft and bibliography file. |

## Useful Starting Code

The most relevant lecture files are:

- `structopt_exercises_lectures/chapter 9/ex1a.m`: MBB topology optimization with OC and sensitivity or density filtering.
- `structopt_exercises_lectures/chapter 9/ex1b.m`: same MBB setup with density and displacement plots.
- `structopt_exercises_lectures/chapter 9/ex3b.m`: MBB topology optimization using GCMMA.
- `structopt_exercises_lectures/chapter 9/ex4.m`: MBB topology optimization using MMA with sensitivity filtering, density filtering, or Heaviside projection.
- `structopt/mma.m`: MMA wrapper around Svanberg's implementation.
- `structopt/gcmma.m`: GCMMA wrapper used in the lecture examples.
- `structopt/plottaylorfdcheck.m`, `structopt/printcomponentfdcheck.m`, `structopt/printdirectionalfdcheck.m`: helpers for sensitivity verification.

The lecture MBB scripts are script-style examples. For the assignment, the new code should be factored into reusable functions so each experiment can be run reproducibly and each figure in the report can be regenerated.

## MATLAB Workspace

The assignment workspace has been started in `matlab/`:

| Path | Purpose |
| --- | --- |
| `matlab/run_assignment.m` | Entry point for the assignment workflow. It sets paths, loads the assignment configuration, lists the planned experiments, and runs implemented steps. |
| `matlab/src/setup_project.m` | Adds local MATLAB paths and creates output folders. |
| `matlab/src/assignment_config.m` | Central configuration for mesh sizes, filters, optimizer settings, and planned experiments. |
| `matlab/src/run_classical_oc_sensitivity.m` | Classical OC sensitivity-filter baseline for assignment step 1. |
| `matlab/src/run_mtop_oc_sensitivity.m` | MTOP OC sensitivity-filter run and comparison against the step 1 baseline. |
| `matlab/figures/` | Target folder for report-ready figures. |
| `matlab/results/` | Target folder for `.mat` result files, timing histories, and intermediate data. |
| `matlab/readme.txt` | Short MATLAB-file overview intended to evolve into the submission ZIP readme. |

From MATLAB, run:

```matlab
cd('C:\Users\campa\Desktop\Structural-Optimization-LUKASTIEBERT\matlab')
run_assignment
```

The current workflow runs assignment steps 1 and 2 and generates the classical OC baseline, MTOP OC result, design difference figure, and convergence comparison figures. The density filtering, MMA, Heaviside projection, and sensitivity verification experiments still need to be implemented in later steps.

## Implementation Plan

Recommended function split for the assignment code:

- `run_classical_oc.m`: classical OC baseline and density-filter variant.
- `run_mtop_oc.m`: MTOP OC runner with coarse analysis mesh and fine density mesh.
- `run_mtop_mma.m`: MTOP MMA runner for sensitivity filtering, density filtering, and Heaviside projection.
- `mbb_fe_model.m`: finite element model construction for the half MBB beam.
- `mtop_mapping.m`: density-to-analysis mapping for `5 x 5` design subcells per finite element.
- `apply_filter.m`: sensitivity filter, density filter, and Heaviside projection.
- `objective_compliance.m`: compliance, volume fraction, and sensitivities.
- `save_result.m` and plotting helpers: save `.mat` histories and generate report figures.
- `check_sensitivities.m`: finite-difference verification for representative design variables.

Key consistency checks:

- Classical baseline and MTOP should both use a `600 x 200` density resolution when designs are visually compared.
- Runtime plots should use wall-clock time from the same machine and the same stopping criteria where possible.
- Store per-iteration compliance, physical volume fraction, change, optimizer status, and elapsed time.
- Use the same color scale and plotting extent for classical and MTOP density plots.
- For MTOP, document how fine-scale densities are interpolated or averaged into analysis-element stiffnesses and how sensitivities are mapped back to density variables.

## Report Checklist

The report should include:

- Problem statement and motivation for MTOP.
- Geometry, supports, loading, material interpolation, and boundary conditions.
- Optimization formulation: design variables, objective, constraints, filters, and projection.
- Classical OC baseline result and convergence history.
- MTOP OC result and comparison against the classical baseline.
- Density-filter comparison.
- MMA results for sensitivity filtering, density filtering, and Heaviside projection thresholds.
- Sensitivity verification using finite differences.
- Runtime comparison in both iteration count and elapsed time.
- Physical interpretation of the final layouts and parameter effects.
- Transparent statement describing GenAI use.
- Bibliography, including Nguyen et al. (2010).

The LaTeX report skeleton is available in `report/main.tex`.
