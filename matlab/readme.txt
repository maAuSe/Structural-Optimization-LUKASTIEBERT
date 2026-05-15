MTOP assignment, MATLAB workspace

How to run:

Open MATLAB in this folder and run run_assignment.m. It sets the paths,
runs the finite-difference sensitivity check, then runs every experiment
end to end. Per-iteration history (.csv), result struct (.mat) and text
summary land in results/; PNG and PDF figures land in figures/.

Cached results in results/*.mat are reused on rerun. Delete a .mat file
to force its step to rebuild.

Files in src: 
setup_project.m         Adds paths, makes the output folders.
assignment_config.m     All shared parameters and the experiment list.

run_classical_oc_sensitivity.m   Step 1. Classical 600 x 200 OC, sens. filter.
run_mtop_oc_sensitivity.m        Step 2. MTOP OC, 120 x 40 FE with 5 x 5 cells.
run_classical_oc_density.m       Step 3a. Classical OC with density filter.
run_mtop_oc_density.m            Step 3b. MTOP OC with density filter.
run_mtop_mma_sensitivity.m       Step 4a. MTOP MMA, sens. filter.
run_mtop_mma_density.m           Step 4b. MTOP MMA, density filter.
run_mtop_mma_heaviside.m         Step 4c. MTOP MMA, density filter + Heaviside.
Called for eta = 0.3, 0.5, 0.7.

mtop_subcell_stiffness.m  The 8 x 8 x 25 sub-cell stiffness templates for
                          the Nguyen Q4/n25 multiresolution assembly.
verify_sensitivities.m    Finite-difference check of dC/drho and dV/drho
                          for the classical, MTOP, and full filter +
                          Heaviside chains. Stored in
                          results/sensitivity_verification.mat.

export_*.m                Standalone report figures (SIMP curve, cone
                          kernel, runtime summary, phase breakdown,
                          density-filter sensitivity fields).

Each runner is self-contained: the FE assembly, filter, OC or MMA update,
summary writer, and figure exporter all live as local functions in the
same file.

Convergence stopping:

The density-filtered runs (classical and MTOP, OC and MMA) do not satisfy
the design-change tolerance epsilon = 0.01 cleanly: the change oscillates
in the 0.04 to 0.07 range because the multiplicative OC update does not
damp the high-frequency content of the filtered design, and MMA inherits
the same behavior. To stop cleanly we add a secondary criterion: the run
terminates when the relative spread of the compliance over a 20-iteration
window falls below 1e-4. With that, the density-filtered OC runs stop at
237, 299, 419 iterations (max cap 500). The MMA + Heaviside runs do reach
epsilon = 0.01 once beta = beta_max = 16, in 419, 563, 522 iterations for
eta = 0.3, 0.5, 0.7 (max cap 800).

Timing:

Every iteration logs wall-clock time per phase (FE assembly, FE solve,
filter / chain rule, optimizer). Totals are written to the summary text
file and to the per-iteration CSV. This is what the report uses to show
where the MTOP speedup comes from (FE solve dominates the gain, the
filter and optimizer scale with the density mesh and partly absorb it).
