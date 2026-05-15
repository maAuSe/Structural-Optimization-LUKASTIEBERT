function cfg = assignment_config(paths)

if nargin < 1 || isempty(paths)
  paths = setup_project();
end

parameters = struct();
parameters.volfrac = 0.5;
parameters.penal = 3;
parameters.classicalNelx = 600;
parameters.classicalNely = 200;
parameters.classicalFilterRadius = 24;
parameters.mtopNelx = 120;
parameters.mtopNely = 40;
parameters.densityPerElement = [5 5];
parameters.mtopDensityNelx = parameters.mtopNelx * parameters.densityPerElement(1);
parameters.mtopDensityNely = parameters.mtopNely * parameters.densityPerElement(2);
parameters.tol = 0.01;
parameters.move = 0.2;
parameters.heavisideBeta0 = 1;
parameters.heavisideBetaMax = 16;
parameters.heavisideEtaValues = [0.3 0.5 0.7];

experiments = struct('id', {}, 'approach', {}, 'optimizer', {}, 'filter', {}, ...
  'feNelx', {}, 'feNely', {}, 'densityNelx', {}, 'densityNely', {}, ...
  'filterRadius', {}, 'eta', {});

experiments(end+1) = make_experiment('classical_oc_sensitivity', 'classical', 'oc', 'sensitivity', ...
  parameters.classicalNelx, parameters.classicalNely, ...
  parameters.classicalNelx, parameters.classicalNely, parameters.classicalFilterRadius, NaN);

experiments(end+1) = make_experiment('mtop_oc_sensitivity', 'mtop', 'oc', 'sensitivity', ...
  parameters.mtopNelx, parameters.mtopNely, ...
  parameters.mtopDensityNelx, parameters.mtopDensityNely, parameters.classicalFilterRadius, NaN);

experiments(end+1) = make_experiment('classical_oc_density', 'classical', 'oc', 'density', ...
  parameters.classicalNelx, parameters.classicalNely, ...
  parameters.classicalNelx, parameters.classicalNely, parameters.classicalFilterRadius, NaN);

experiments(end+1) = make_experiment('mtop_oc_density', 'mtop', 'oc', 'density', ...
  parameters.mtopNelx, parameters.mtopNely, ...
  parameters.mtopDensityNelx, parameters.mtopDensityNely, parameters.classicalFilterRadius, NaN);

experiments(end+1) = make_experiment('mtop_mma_sensitivity', 'mtop', 'mma', 'sensitivity', ...
  parameters.mtopNelx, parameters.mtopNely, ...
  parameters.mtopDensityNelx, parameters.mtopDensityNely, parameters.classicalFilterRadius, NaN);

experiments(end+1) = make_experiment('mtop_mma_density', 'mtop', 'mma', 'density', ...
  parameters.mtopNelx, parameters.mtopNely, ...
  parameters.mtopDensityNelx, parameters.mtopDensityNely, parameters.classicalFilterRadius, NaN);

for k = 1:numel(parameters.heavisideEtaValues)
  eta = parameters.heavisideEtaValues(k);
  id = sprintf('mtop_mma_heaviside_eta_%03i', round(100 * eta));
  experiments(end+1) = make_experiment(id, 'mtop', 'mma', 'heaviside', ...
    parameters.mtopNelx, parameters.mtopNely, ...
    parameters.mtopDensityNelx, parameters.mtopDensityNely, parameters.classicalFilterRadius, eta);
end

cfg = struct();
cfg.paths = paths;
cfg.parameters = parameters;
cfg.experiments = experiments;

end

function experiment = make_experiment(id, approach, optimizer, filter, feNelx, feNely, densityNelx, densityNely, filterRadius, eta)
experiment = struct();
experiment.id = id;
experiment.approach = approach;
experiment.optimizer = optimizer;
experiment.filter = filter;
experiment.feNelx = feNelx;
experiment.feNely = feNely;
experiment.densityNelx = densityNelx;
experiment.densityNely = densityNely;
experiment.filterRadius = filterRadius;
experiment.eta = eta;
end
