MATLAB files for the Structural Optimization MTOP assignment

Entry point:
- run_assignment.m
  Sets up paths, loads the assignment configuration, runs the sensitivity
  verification, and runs the implemented assignment steps. Existing step 1
  results are reused when available; step 2 is always executed against
  the freshly assembled MTOP stiffness.

Source folder (matlab/src):
- setup_project.m
  Adds local helper paths and ensures the figures and results folders exist.

- assignment_config.m
  Stores the main assignment parameters and planned experiment cases:
  classical OC baseline, MTOP OC, density filtering, MMA, and Heaviside
  projection threshold variants.

- run_classical_oc_sensitivity.m
  Reproducible runner for assignment step 1. Translates the lecture MBB
  beam OC code (chapter 9 / ex1a.m) into a runner with logging, history
  storage and figure export, for a 600 x 200 mesh, volume fraction 0.5,
  SIMP penalization 3, and sensitivity filtering with radius 24 elements.

- run_mtop_oc_sensitivity.m
  Reproducible runner for assignment step 2. Solves the same MBB problem
  with a 120 x 40 finite element mesh and a 600 x 200 density mesh, using
  5 x 5 density cells per finite element. The element stiffness matrix is
  assembled by the Q4/n25 multiresolution scheme of Nguyen et al. (2010),
  Eq. (11): each density cell contributes a precomputed B^T D^0 B template
  weighted by the SIMP-penalized cell density.

- run_classical_oc_density.m
  Reproducible runner for the classical half of assignment step 3.
  Same problem and parameters as run_classical_oc_sensitivity, but with
  a density filter (the design variables are convolved with the cone
  kernel before SIMP analysis, and the chain rule is applied to the
  sensitivities).

- run_mtop_oc_density.m
  Reproducible runner for the MTOP half of assignment step 3.
  Same MTOP setup as run_mtop_oc_sensitivity, but with a density filter
  applied on the 600 x 200 density mesh.

- run_mtop_mma_sensitivity.m
- run_mtop_mma_density.m
  Reproducible runners for the MMA-based variants of step 4. Each runner
  replaces the OC update by a call to the structopt MMA implementation
  (mma.m) while keeping the multiresolution assembly and the corresponding
  filter chain rule unchanged.

- run_mtop_mma_heaviside.m
  Reproducible runner for the MMA + Heaviside-projection variant of
  step 4. Takes the projection threshold eta as a second argument; the
  Heaviside continuation parameter beta is doubled from 1 to 16 every
  50 iterations following Wang, Lazarov and Sigmund (2011).
  run_assignment.m calls this runner three times with eta in {0.3, 0.5, 0.7}.

- mtop_subcell_stiffness.m
  Helper that returns the 8 x 8 x 25 array of sub-cell stiffness templates
  used by the MTOP assembly.

- verify_sensitivities.m
  Finite-difference verification of the compliance and volume sensitivities
  for the classical and MTOP formulations, plus a check of the full filter
  + Heaviside + multiresolution chain rule used in step 4. The maximum
  absolute and relative errors are reported and stored in
  results/sensitivity_verification.mat. The verification is invoked by
  run_assignment.m before any production runs.

Output folders:
- figures/
  Generated report figures.
  Step 1: classical_oc_sensitivity_design.{png,pdf},
          classical_oc_sensitivity_convergence_{iter,time}.{png,pdf}.
  Step 2: mtop_oc_sensitivity_design.{png,pdf},
          mtop_oc_sensitivity_design_difference.{png,pdf},
          mtop_oc_sensitivity_convergence_{iter,time}_compare.{png,pdf}.
  Step 3: classical_oc_density_design.{png,pdf},
          classical_oc_density_convergence_{iter,time}.{png,pdf},
          mtop_oc_density_design.{png,pdf},
          mtop_oc_density_design_difference.{png,pdf},
          mtop_oc_density_convergence_{iter,time}_compare.{png,pdf}.
  Step 4: mtop_mma_sensitivity_design.{png,pdf},
          mtop_mma_sensitivity_convergence_{iter,time}.{png,pdf},
          mtop_mma_density_design.{png,pdf},
          mtop_mma_density_convergence_{iter,time}.{png,pdf},
          mtop_mma_heaviside_eta_{030,050,070}_design.{png,pdf},
          mtop_mma_heaviside_eta_{030,050,070}_convergence_{iter,time}.{png,pdf}.

- results/
  Generated result files for each runner: a .mat file with the full
  result struct, a .csv file with the per-iteration history (iteration
  counter, elapsed wall-clock time, compliance, volume fraction,
  maximum design change, and beta for the Heaviside runs) and a .txt
  summary with the key scalars. The verification report is stored in
  sensitivity_verification.mat.

Implementation status:
All four assignment steps are implemented, including the proper
Nguyen Q4/n25 multiresolution stiffness, a finite-difference
sensitivity verification routine for the classical, MTOP and full
filter + Heaviside chains, and the MMA + Heaviside continuation
strategy of Wang, Lazarov and Sigmund (2011).

Notes on convergence:
The density-filtered runs (classical OC and MTOP OC, step 3, and MTOP
MMA density, step 4) do not satisfy the design-change tolerance
epsilon = 0.01: the maximum design change oscillates between 0.04 and
0.07 because the multiplicative OC update does not damp the
high-frequency components of the update under density filtering, and
MMA inherits the same pathology. To stop these runs cleanly we use a
secondary compliance-stability criterion: the optimization terminates
when the relative spread of the compliance over a 20-iteration window
drops below 1e-4. With this criterion all three density-filtered runs
stop at well-defined iteration counts (237, 299, 419) rather than at
the maximum-iteration cap (now set to 500). The MMA + Heaviside runs
of step 4 produce essentially black-and-white designs under continuation
on beta and all three eta values (0.3, 0.5, 0.7) satisfy the
design-change tolerance once beta = beta_max = 16, in 419, 563 and
522 iterations respectively (max iter cap = 800 for the Heaviside runs).

Notes on timing:
Each runner is instrumented with per-phase wall-clock timing for the
FE assembly, FE solve, filter / chain rule, and optimizer (OC bisection
or MMA subproblem). The cumulative phase times are written to the
summary file and to the per-iteration history CSV, in addition to the
total elapsed time. This is used in the report to decompose the
measured speedups (FE solve dominates the gain; the filter and the
optimizer scale with the density mesh and absorb part of it).

Reproducing the figures:
- Open MATLAB in matlab/ and run run_assignment.m. Step 1 and the
  classical half of step 3 are skipped if their .mat files exist in
  results/. To rerun a step from scratch, delete its .mat file
  beforehand. All other steps always run.
