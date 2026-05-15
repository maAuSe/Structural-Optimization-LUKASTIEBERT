% entry pointttt for the workspace of the MASTERCLASS made by Tiebert & Lukas
%
% this script sets up paths, creates output folders, prints the experiment
% matrix encoded in src/assignment_config.m, and runs the implemented steps like a master.

clear;
clc;

matlabDir = fileparts(mfilename('fullpath'));
addpath(fullfile(matlabDir, 'src'));

paths = setup_project();
cfg = assignment_config(paths);

fprintf('\nStructural Optimization MTOP assignment\n');
fprintf('Repository root: %s\n', cfg.paths.root);
fprintf('Results folder:  %s\n', cfg.paths.results);
fprintf('Figures folder:  %s\n\n', cfg.paths.figures);

fprintf('Core parameters\n');
fprintf('  Volume fraction: %.3f\n', cfg.parameters.volfrac);
fprintf('  Penalization:    %.3f\n', cfg.parameters.penal);
fprintf('  Classical mesh:  %i x %i\n', cfg.parameters.classicalNelx, cfg.parameters.classicalNely);
fprintf('  MTOP FE mesh:    %i x %i\n', cfg.parameters.mtopNelx, cfg.parameters.mtopNely);
fprintf('  MTOP density:    %i x %i per FE\n\n', cfg.parameters.densityPerElement(1), cfg.parameters.densityPerElement(2));

fprintf('Planned experiments\n');
for k = 1:numel(cfg.experiments)
  e = cfg.experiments(k);
  fprintf('%2i. %-30s | %-9s | %-3s | %-11s | FE %4i x %-4i | density %4i x %-4i\n', ...
    k, e.id, e.approach, e.optimizer, e.filter, e.feNelx, e.feNely, e.densityNelx, e.densityNely);
end

fprintf('\nExecuting implemented assignment steps.\n');

verifyReport = verify_sensitivities(cfg);

classicalPath = fullfile(cfg.paths.results, 'classical_oc_sensitivity.mat');
if isfile(classicalPath)
  data = load(classicalPath, 'result');
  classicalResult = data.result;
  fprintf('Step 1 result found and reused: %s\n', classicalPath);
else
  classicalResult = run_classical_oc_sensitivity(cfg);
end

fprintf('Step 1 completed.\n');
fprintf('  Iterations:       %i\n', classicalResult.iterations);
fprintf('  Elapsed time [s]: %.2f\n', classicalResult.elapsedSeconds);
fprintf('  Compliance:       %.6f\n', classicalResult.finalCompliance);
fprintf('  Volume fraction:  %.6f\n', classicalResult.finalVolumeFraction);

mtopResult = run_mtop_oc_sensitivity(cfg, classicalResult);

fprintf('Step 2 completed.\n');
fprintf('  Iterations:              %i\n', mtopResult.iterations);
fprintf('  Elapsed time [s]:        %.2f\n', mtopResult.elapsedSeconds);
fprintf('  MTOP compliance:         %.6f\n', mtopResult.finalCompliance);
fprintf('  Fine-mesh compliance:    %.6f\n', mtopResult.fineMeshCompliance);
fprintf('  Volume fraction:         %.6f\n', mtopResult.finalVolumeFraction);
fprintf('  Speedup vs classical:    %.2f\n', mtopResult.speedupAgainstClassical);

coarsePath = fullfile(cfg.paths.results, 'classical_oc_coarse_sensitivity.mat');
if isfile(coarsePath)
  data = load(coarsePath, 'result');
  coarseResult = data.result;
  fprintf('Coarse-coarse sensitivity result found and reused: %s\n', coarsePath);
else
  coarseResult = run_classical_oc_coarse_sensitivity(cfg);
end

fprintf('Coarse-coarse sensitivity baseline completed.\n');
fprintf('  Iterations:              %i\n', coarseResult.iterations);
fprintf('  Elapsed time [s]:        %.2f\n', coarseResult.elapsedSeconds);
fprintf('  Native coarse compliance: %.6f\n', coarseResult.finalCompliance);
fprintf('  Fine-mesh compliance:    %.6f\n', coarseResult.fineMeshCompliance);
fprintf('  Grayness index:          %.6f\n', coarseResult.grayness);

classicalDensityPath = fullfile(cfg.paths.results, 'classical_oc_density.mat');
if isfile(classicalDensityPath)
  data = load(classicalDensityPath, 'result');
  classicalDensityResult = data.result;
  fprintf('Step 3 (classical density) result found and reused: %s\n', classicalDensityPath);
else
  classicalDensityResult = run_classical_oc_density(cfg);
end

fprintf('Step 3 (classical density) completed.\n');
fprintf('  Iterations:       %i\n', classicalDensityResult.iterations);
fprintf('  Elapsed time [s]: %.2f\n', classicalDensityResult.elapsedSeconds);
fprintf('  Compliance:       %.6f\n', classicalDensityResult.finalCompliance);
fprintf('  Volume fraction:  %.6f\n', classicalDensityResult.finalVolumeFraction);

mtopDensityResult = run_mtop_oc_density(cfg, classicalDensityResult);

fprintf('Step 3 (MTOP density) completed.\n');
fprintf('  Iterations:              %i\n', mtopDensityResult.iterations);
fprintf('  Elapsed time [s]:        %.2f\n', mtopDensityResult.elapsedSeconds);
fprintf('  MTOP compliance:         %.6f\n', mtopDensityResult.finalCompliance);
fprintf('  Fine-mesh compliance:    %.6f\n', mtopDensityResult.fineMeshCompliance);
fprintf('  Volume fraction:         %.6f\n', mtopDensityResult.finalVolumeFraction);
fprintf('  Speedup vs classical:    %.2f\n', mtopDensityResult.speedupAgainstClassical);

mmaSensitivityResult = run_mtop_mma_sensitivity(cfg);

fprintf('Step 4 (MTOP MMA sensitivity) completed.\n');
fprintf('  Iterations:           %i\n', mmaSensitivityResult.iterations);
fprintf('  Elapsed time [s]:     %.2f\n', mmaSensitivityResult.elapsedSeconds);
fprintf('  MTOP compliance:      %.6f\n', mmaSensitivityResult.finalCompliance);
fprintf('  Fine-mesh compliance: %.6f\n', mmaSensitivityResult.fineMeshCompliance);
fprintf('  Volume fraction:      %.6f\n', mmaSensitivityResult.finalVolumeFraction);

mmaDensityResult = run_mtop_mma_density(cfg);

fprintf('Step 4 (MTOP MMA density) completed.\n');
fprintf('  Iterations:           %i\n', mmaDensityResult.iterations);
fprintf('  Elapsed time [s]:     %.2f\n', mmaDensityResult.elapsedSeconds);
fprintf('  MTOP compliance:      %.6f\n', mmaDensityResult.finalCompliance);
fprintf('  Fine-mesh compliance: %.6f\n', mmaDensityResult.fineMeshCompliance);
fprintf('  Volume fraction:      %.6f\n', mmaDensityResult.finalVolumeFraction);

heavisideEtaValues = cfg.parameters.heavisideEtaValues;
heavisideResults = cell(1, numel(heavisideEtaValues));
for k = 1:numel(heavisideEtaValues)
  eta = heavisideEtaValues(k);
  heavisideResults{k} = run_mtop_mma_heaviside(cfg, eta);
  fprintf('Step 4 (MTOP MMA Heaviside, eta = %.2f) completed.\n', eta);
  fprintf('  Iterations:           %i\n', heavisideResults{k}.iterations);
  fprintf('  Elapsed time [s]:     %.2f\n', heavisideResults{k}.elapsedSeconds);
  fprintf('  Final beta:           %.4f\n', heavisideResults{k}.finalBeta);
  fprintf('  MTOP compliance:      %.6f\n', heavisideResults{k}.finalCompliance);
  fprintf('  Fine-mesh compliance: %.6f\n', heavisideResults{k}.fineMeshCompliance);
  fprintf('  Volume fraction:      %.6f\n', heavisideResults{k}.finalVolumeFraction);
end

fprintf('Sensitivity verification:\n');
fprintf('  Classical (dC/drho)        max rel FD error: %.3e\n', verifyReport.classical.dcMaxRelError);
fprintf('  MTOP      (dC/drho)        max rel FD error: %.3e\n', verifyReport.mtop.dcMaxRelError);
fprintf('  Heaviside (dC/dx full chain) max rel FD error: %.3e\n', verifyReport.heaviside.dcMaxRelError);

export_runtime_summary_chart();
export_phase_breakdown_chart();
export_density_sensitivity_figures(cfg);
