%RUN_ASSIGNMENT Entry point for the MTOP assignment workspace.
%
% This script sets up paths, creates output folders, prints the experiment
% matrix encoded in src/assignment_config.m, and runs the implemented steps.

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
