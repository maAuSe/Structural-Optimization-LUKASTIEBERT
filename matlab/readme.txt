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

- mtop_subcell_stiffness.m
  Helper that returns the 8 x 8 x 25 array of sub-cell stiffness templates
  used by the MTOP assembly.

- verify_sensitivities.m
  Finite-difference verification of the compliance and volume sensitivities
  of both formulations on a small MBB instance. Reports the maximum
  absolute and relative errors and stores them in
  results/sensitivity_verification.mat. The verification is invoked by
  run_assignment.m before step 1 is started.

Output folders:
- figures/
  Generated report figures. Step 1 and step 2 outputs:
  classical_oc_sensitivity_design.png/.pdf
  classical_oc_sensitivity_convergence_iter.png/.pdf
  classical_oc_sensitivity_convergence_time.png/.pdf
  mtop_oc_sensitivity_design.png/.pdf
  mtop_oc_sensitivity_design_difference.png/.pdf
  mtop_oc_sensitivity_convergence_iter_compare.png/.pdf
  mtop_oc_sensitivity_convergence_time_compare.png/.pdf

- results/
  Generated result files. Step 1 and step 2 outputs:
  classical_oc_sensitivity.mat
  classical_oc_sensitivity_history.csv
  classical_oc_sensitivity_summary.txt
  mtop_oc_sensitivity.mat
  mtop_oc_sensitivity_history.csv
  mtop_oc_sensitivity_summary.txt
  sensitivity_verification.mat

Implementation status:
Assignment steps 1 and 2 are implemented, including the proper Nguyen
Q4/n25 multiresolution stiffness and a finite-difference sensitivity
verification routine. The remaining density filtering, MMA and Heaviside
projection experiments still need to be implemented.

Reproducing the figures:
- Open MATLAB in matlab/ and run run_assignment.m. Step 1 is skipped if
  classical_oc_sensitivity.mat exists in results/. To rerun step 1 from
  scratch, delete that file beforehand. Step 2 always runs because its
  cache was invalidated by the move from the averaged-density assembly to
  the proper Nguyen Q4/n25 assembly.
