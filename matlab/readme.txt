MATLAB files for the Structural Optimization MTOP assignment

Entry point:
- run_assignment.m
  Sets up paths, loads the assignment configuration, creates output folders,
  prints the planned experiment matrix, and runs the implemented assignment
  steps. Existing step 1 results are reused when available.

Source folder:
- src/setup_project.m
  Adds local helper paths and ensures the figures and results folders exist.

- src/assignment_config.m
  Stores the main assignment parameters and planned experiment cases:
  classical OC baseline, MTOP OC, density filtering, MMA, and Heaviside
  projection threshold variants.

- src/run_classical_oc_sensitivity.m
  Reproducible runner for assignment step 1. It adapts the lecture MBB beam
  OC code for a 600 x 200 mesh, volume fraction 0.5, SIMP penalization 3,
  and sensitivity filtering with radius 24 elements. It stores convergence
  histories, timing data, and report figures.

- src/run_mtop_oc_sensitivity.m
  Reproducible runner for assignment step 2. It solves the same MBB problem
  with a 120 x 40 finite element mesh and a 600 x 200 density mesh, using
  5 x 5 density cells per finite element. It stores the MTOP history, final
  design, comparison metrics against step 1, and report figures.

Output folders:
- figures/
  Generated report figures. Current step 1 and step 2 outputs:
  classical_oc_sensitivity_design.png/.pdf
  classical_oc_sensitivity_convergence_iter.png/.pdf
  classical_oc_sensitivity_convergence_time.png/.pdf
  mtop_oc_sensitivity_design.png/.pdf
  mtop_oc_sensitivity_design_difference.png/.pdf
  mtop_oc_sensitivity_convergence_iter_compare.png/.pdf
  mtop_oc_sensitivity_convergence_time_compare.png/.pdf

- results/
  Generated result files. Current step 1 and step 2 outputs:
  classical_oc_sensitivity.mat
  classical_oc_sensitivity_history.csv
  classical_oc_sensitivity_summary.txt
  step1_run.log
  mtop_oc_sensitivity.mat
  mtop_oc_sensitivity_history.csv
  mtop_oc_sensitivity_summary.txt
  step2_run.log

Implementation status:
Assignment steps 1 and 2 are implemented and were run successfully. The
remaining density filtering, MMA, Heaviside projection, and sensitivity
verification experiments still need to be implemented.
