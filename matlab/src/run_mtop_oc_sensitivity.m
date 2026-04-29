function result = run_mtop_oc_sensitivity(cfg, classicalResult)
%RUN_MTOP_OC_SENSITIVITY MTOP OC sensitivity-filter study for assignment step 2.
%
% The density mesh has 5 x 5 density cells per finite element. The element
% stiffness matrix is assembled following the multiresolution scheme of
% Nguyen et al. (2010), Eq. (11): each density cell contributes a
% precomputed template I_k (B^T D^0 B evaluated at the cell centre, scaled
% by the cell area), weighted by its SIMP-penalised density. The compliance
% sensitivity for a density cell is therefore the local strain-energy
% density of that cell rather than the element-averaged value.

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

if nargin < 2 || isempty(classicalResult)
  classicalResult = load_or_run_classical_reference(cfg);
end

p = cfg.parameters;

problem = struct();
problem.id = 'mtop_oc_sensitivity';
problem.feNelx = p.mtopNelx;
problem.feNely = p.mtopNely;
problem.densityPerElement = p.densityPerElement;
problem.densityNelx = p.mtopDensityNelx;
problem.densityNely = p.mtopDensityNely;
problem.volfrac = p.volfrac;
problem.penal = p.penal;
problem.rmin = p.classicalFilterRadius;
problem.tol = p.tol;
problem.move = p.move;
problem.maxIter = 300;
problem.E0 = 1;
problem.Emin = 1e-9;
problem.nu = 0.3;

fprintf('\nRunning assignment step 2: MTOP OC sensitivity baseline\n');
fprintf('  FE mesh:       %i x %i elements\n', problem.feNelx, problem.feNely);
fprintf('  Density mesh:  %i x %i cells\n', problem.densityNelx, problem.densityNely);
fprintf('  Cells per FE:  %i x %i\n', problem.densityPerElement(1), problem.densityPerElement(2));
fprintf('  Volume target: %.4f\n', problem.volfrac);
fprintf('  Penalization:  %.4f\n', problem.penal);
fprintf('  Filter radius: %.4f density cells\n\n', problem.rmin);

result = solve_mtop_mbb_oc(problem);
result.problem = problem;

densityDifference = result.xPhys - classicalResult.xPhys;
result.reference = make_reference_summary(classicalResult);
result.difference = struct();
result.difference.meanAbsDensity = mean(abs(densityDifference(:)));
result.difference.rmsDensity = sqrt(mean(densityDifference(:).^2));
result.difference.maxAbsDensity = max(abs(densityDifference(:)));
result.fineMeshCompliance = evaluate_fine_mbb_compliance(result.xPhys, problem);
result.fineMeshComplianceRatio = result.fineMeshCompliance / classicalResult.finalCompliance;
result.speedupAgainstClassical = classicalResult.elapsedSeconds / result.elapsedSeconds;

resultPath = fullfile(cfg.paths.results, [problem.id '.mat']);
historyPath = fullfile(cfg.paths.results, [problem.id '_history.csv']);
summaryPath = fullfile(cfg.paths.results, [problem.id '_summary.txt']);

save(resultPath, 'result', '-v7.3');
write_history_csv(historyPath, result.history);
write_summary(summaryPath, result);
export_mtop_oc_figures(cfg.paths.figures, result, classicalResult, densityDifference);

fprintf('\nSaved step 2 result data:\n');
fprintf('  %s\n', resultPath);
fprintf('  %s\n', historyPath);
fprintf('  %s\n', summaryPath);
fprintf('Saved step 2 figures in:\n');
fprintf('  %s\n\n', cfg.paths.figures);

end

function classicalResult = load_or_run_classical_reference(cfg)

resultPath = fullfile(cfg.paths.results, 'classical_oc_sensitivity.mat');
if isfile(resultPath)
  data = load(resultPath, 'result');
  classicalResult = data.result;
else
  fprintf('Classical reference result not found. Running assignment step 1 first.\n');
  classicalResult = run_classical_oc_sensitivity(cfg);
end

end

function result = solve_mtop_mbb_oc(problem)

feNelx = problem.feNelx;
feNely = problem.feNely;
densityPerX = problem.densityPerElement(1);
densityPerY = problem.densityPerElement(2);
densityNelx = problem.densityNelx;
densityNely = problem.densityNely;
volfrac = problem.volfrac;
penal = problem.penal;
rmin = problem.rmin;
tol = problem.tol;
move = problem.move;
maxIter = problem.maxIter;
E0 = problem.E0;
Emin = problem.Emin;
nu = problem.nu;

I_cells = mtop_subcell_stiffness(densityPerX, densityPerY, nu);
nodenrs = reshape(1:(1 + feNelx) * (1 + feNely), 1 + feNely, 1 + feNelx);
edofVec = reshape(2 * nodenrs(1:end-1, 1:end-1) + 1, feNelx * feNely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * feNely + [2 3 0 1] -2 -1], feNelx * feNely, 1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * feNelx * feNely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * feNelx * feNely, 1);

F = sparse(2, 1, -1, 2 * (feNely + 1) * (feNelx + 1), 1);
U = zeros(2 * (feNely + 1) * (feNelx + 1), 1);
fixeddofs = union(1:2:2 * (feNely + 1), 2 * (feNelx + 1) * (feNely + 1));
alldofs = 1:2 * (feNely + 1) * (feNelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

[dy, dx] = meshgrid(-ceil(rmin) + 1:ceil(rmin) - 1, -ceil(rmin) + 1:ceil(rmin) - 1);
h = max(0, rmin - sqrt(dx.^2 + dy.^2));
Hs = conv2(ones(densityNely, densityNelx), h, 'same');

x = repmat(volfrac, densityNely, densityNelx);
xPhys = x;

history = struct();
history.iteration = zeros(maxIter, 1);
history.timeSeconds = zeros(maxIter, 1);
history.compliance = zeros(maxIter, 1);
history.volumeFraction = zeros(maxIter, 1);
history.change = zeros(maxIter, 1);

timer = tic;
change = 1;
iter = 0;

while change > tol && iter < maxIter
  iter = iter + 1;

  xPhys = x;
  sK = mtop_assemble_stiffness(xPhys, I_cells, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);

  K = sparse(iK, jK, sK);
  K = (K + K') / 2;
  U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

  ceCells = mtop_strain_energy(U(edofMat), I_cells, ...
    densityPerX, densityPerY, feNelx, feNely);
  c = sum(sum((Emin + (E0 - Emin) * xPhys.^penal) .* ceCells));
  v = mean(xPhys(:));

  dcdx = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ceCells;
  dvdx = ones(densityNely, densityNelx) / densityNelx / densityNely;
  dcdx = conv2(dcdx .* xPhys, h, 'same') ./ Hs ./ max(1e-3, xPhys);

  l1 = 0;
  l2 = 1e9;
  while (l2 - l1) / (l1 + l2) > 1e-3
    lmid = 0.5 * (l2 + l1);
    xCandidate = x .* sqrt(-dcdx ./ dvdx / lmid);
    xnew = max(0, max(x - move, min(1, min(x + move, xCandidate))));
    if sum(xnew(:)) > volfrac * densityNelx * densityNely
      l1 = lmid;
    else
      l2 = lmid;
    end
  end

  change = max(abs(xnew(:) - x(:)));
  elapsed = toc(timer);

  history.iteration(iter) = iter;
  history.timeSeconds(iter) = elapsed;
  history.compliance(iter) = c;
  history.volumeFraction(iter) = v;
  history.change(iter) = change;

  fprintf('[%s] iter: %4i | change: %9.5f | c: %12.5f | v: %9.5f | time: %9.2f s\n', ...
    datestr(now, 'HH:MM:SS'), iter, change, c, v, elapsed);

  if change > tol
    x = xnew;
  end
end

fields = fieldnames(history);
for k = 1:numel(fields)
  history.(fields{k}) = history.(fields{k})(1:iter);
end

result = struct();
result.id = problem.id;
result.xDesign = x;
result.xPhys = xPhys;
result.U = U;
result.history = history;
result.finalCompliance = history.compliance(end);
result.finalVolumeFraction = history.volumeFraction(end);
result.finalChange = history.change(end);
result.iterations = iter;
result.elapsedSeconds = history.timeSeconds(end);
result.converged = change <= tol;
result.filterKernel = h;
result.filterWeights = Hs;

if ~result.converged
  warning('run_mtop_oc_sensitivity:notConverged', ...
    'Maximum iteration count reached before the change tolerance was met.');
end

end

function sK = mtop_assemble_stiffness(xPhys, I_cells, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin)
%MTOP_ASSEMBLE_STIFFNESS Element stiffness entries for the Q4/nN scheme.
%
% Returns a (64 * feNelx * feNely) x 1 vector compatible with the standard
% (iK, jK) sparse-assembly pattern. For each analysis element e and each
% sub-cell k inside it,
%
%   K_e = sum_k (E_min + (E_0 - E_min) * rho_k^p) * I_cells(:, :, k).

nFe = feNelx * feNely;
rho4 = reshape(xPhys, densityPerY, feNely, densityPerX, feNelx);
sK = zeros(64 * nFe, 1);

idx = 0;
for sub_x = 1:densityPerX
  for sub_y = 1:densityPerY
    idx = idx + 1;
    rhoSub = reshape(rho4(sub_y, :, sub_x, :), feNely, feNelx);
    Ek = Emin + (E0 - Emin) * rhoSub.^penal;
    Ik = I_cells(:, :, idx);
    sK = sK + reshape(Ik(:) * Ek(:)', 64 * nFe, 1);
  end
end

end


function ceCells = mtop_strain_energy(ueK, I_cells, ...
    densityPerX, densityPerY, feNelx, feNely)
%MTOP_STRAIN_ENERGY Per-density-cell strain-energy density u_e^T I_k u_e.
%
% Returns a (densityPerY * feNely) x (densityPerX * feNelx) matrix giving,
% for every density cell, the value u_e^T I_k u_e where e is the analysis
% element containing the cell and k is the cell's local index inside e.

ceCells = zeros(densityPerY * feNely, densityPerX * feNelx);

idx = 0;
for sub_x = 1:densityPerX
  for sub_y = 1:densityPerY
    idx = idx + 1;
    Ik = I_cells(:, :, idx);
    uIuVec = sum((ueK * Ik) .* ueK, 2);
    uIuGrid = reshape(uIuVec, feNely, feNelx);
    ceCells(sub_y:densityPerY:end, sub_x:densityPerX:end) = uIuGrid;
  end
end

end

function reference = make_reference_summary(classicalResult)

reference = struct();
reference.id = classicalResult.id;
reference.iterations = classicalResult.iterations;
reference.elapsedSeconds = classicalResult.elapsedSeconds;
reference.finalCompliance = classicalResult.finalCompliance;
reference.finalVolumeFraction = classicalResult.finalVolumeFraction;
reference.finalChange = classicalResult.finalChange;

end

function compliance = evaluate_fine_mbb_compliance(xPhys, problem)

nelx = size(xPhys, 2);
nely = size(xPhys, 1);
E0 = problem.E0;
Emin = problem.Emin;
penal = problem.penal;
nu = problem.nu;

KE = element_stiffness_matrix(nu);
nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(2 * nodenrs(1:end-1, 1:end-1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1);

F = sparse(2, 1, -1, 2 * (nely + 1) * (nelx + 1), 1);
U = zeros(2 * (nely + 1) * (nelx + 1), 1);
fixeddofs = union(1:2:2 * (nely + 1), 2 * (nelx + 1) * (nely + 1));
alldofs = 1:2 * (nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

sK = reshape(KE(:) * (Emin + xPhys(:)'.^penal * (E0 - Emin)), 64 * nelx * nely, 1);
K = sparse(iK, jK, sK);
K = (K + K') / 2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
compliance = sum(sum((Emin + xPhys.^penal * (E0 - Emin)) .* ce));

end

function KE = element_stiffness_matrix(nu)

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1 / (1 - nu^2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);

end

function write_history_csv(historyPath, history)

T = table(history.iteration, history.timeSeconds, history.compliance, ...
  history.volumeFraction, history.change, ...
  'VariableNames', {'iteration', 'time_seconds', 'compliance', 'volume_fraction', 'change'});
writetable(T, historyPath);

end

function write_summary(summaryPath, result)

fid = fopen(summaryPath, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'MTOP OC sensitivity baseline\n');
fprintf(fid, 'FE mesh: %i x %i elements\n', result.problem.feNelx, result.problem.feNely);
fprintf(fid, 'Density mesh: %i x %i cells\n', result.problem.densityNelx, result.problem.densityNely);
fprintf(fid, 'Density cells per FE: %i x %i\n', result.problem.densityPerElement(1), result.problem.densityPerElement(2));
fprintf(fid, 'Volume fraction target: %.6f\n', result.problem.volfrac);
fprintf(fid, 'Penalization power: %.6f\n', result.problem.penal);
fprintf(fid, 'Sensitivity filter radius: %.6f density cells\n', result.problem.rmin);
fprintf(fid, 'Move limit: %.6f\n', result.problem.move);
fprintf(fid, 'Change tolerance: %.6f\n', result.problem.tol);
fprintf(fid, 'Converged: %i\n', result.converged);
fprintf(fid, 'Iterations: %i\n', result.iterations);
fprintf(fid, 'Elapsed time [s]: %.6f\n', result.elapsedSeconds);
fprintf(fid, 'Final MTOP compliance: %.12f\n', result.finalCompliance);
fprintf(fid, 'Final fine-mesh compliance of MTOP design: %.12f\n', result.fineMeshCompliance);
fprintf(fid, 'Fine-mesh compliance ratio against classical: %.12f\n', result.fineMeshComplianceRatio);
fprintf(fid, 'Final volume fraction: %.12f\n', result.finalVolumeFraction);
fprintf(fid, 'Final maximum design change: %.12f\n', result.finalChange);
fprintf(fid, 'Speedup against classical wall-clock time: %.12f\n', result.speedupAgainstClassical);
fprintf(fid, 'Mean absolute density difference: %.12f\n', result.difference.meanAbsDensity);
fprintf(fid, 'RMS density difference: %.12f\n', result.difference.rmsDensity);
fprintf(fid, 'Maximum absolute density difference: %.12f\n', result.difference.maxAbsDensity);

end

function export_mtop_oc_figures(figuresDir, result, classicalResult, densityDifference)

if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

history = result.history;
classicalHistory = classicalResult.history;

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(1 - result.xPhys);
colormap(fig, gray);
caxis([0 1]);
axis equal tight off;
export_figure(fig, figuresDir, 'mtop_oc_sensitivity_design');
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
imagesc(1 - classicalResult.xPhys);
colormap(gca, gray);
caxis([0 1]);
axis equal tight off;
title('Classical design');
nexttile;
imagesc(1 - result.xPhys);
colormap(gca, gray);
caxis([0 1]);
axis equal tight off;
title('MTOP design');
nexttile;
imagesc(densityDifference);
colormap(gca, blue_white_red(256));
diffLimit = max(0.05, max(abs(densityDifference(:))));
caxis([-diffLimit diffLimit]);
axis equal tight off;
title('MTOP density minus classical density');
cb = colorbar;
cb.Label.String = 'Density difference';
export_figure(fig, figuresDir, 'mtop_oc_sensitivity_design_difference');
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(classicalHistory.iteration, classicalHistory.compliance, 'LineWidth', 1.4);
hold on;
plot(history.iteration, history.compliance, 'LineWidth', 1.4);
grid on;
xlabel('Iteration number');
ylabel('Compliance');
title('Compliance convergence by iteration');
legend({'Classical', 'MTOP'}, 'Location', 'northeast');
nexttile;
plot(classicalHistory.iteration, classicalHistory.volumeFraction, 'LineWidth', 1.4);
hold on;
plot(history.iteration, history.volumeFraction, 'LineWidth', 1.4);
yline(result.problem.volfrac, ':', 'Target', 'LineWidth', 1.0);
grid on;
xlabel('Iteration number');
ylabel('Volume fraction');
title('Volume convergence by iteration');
legend({'Classical', 'MTOP', 'Target'}, 'Location', 'best');
export_figure(fig, figuresDir, 'mtop_oc_sensitivity_convergence_iter_compare');
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(classicalHistory.timeSeconds, classicalHistory.compliance, 'LineWidth', 1.4);
hold on;
plot(history.timeSeconds, history.compliance, 'LineWidth', 1.4);
grid on;
xlabel('Computation time [s]');
ylabel('Compliance');
title('Compliance convergence by time');
legend({'Classical', 'MTOP'}, 'Location', 'northeast');
nexttile;
plot(classicalHistory.timeSeconds, classicalHistory.volumeFraction, 'LineWidth', 1.4);
hold on;
plot(history.timeSeconds, history.volumeFraction, 'LineWidth', 1.4);
yline(result.problem.volfrac, ':', 'Target', 'LineWidth', 1.0);
grid on;
xlabel('Computation time [s]');
ylabel('Volume fraction');
title('Volume convergence by time');
legend({'Classical', 'MTOP', 'Target'}, 'Location', 'best');
export_figure(fig, figuresDir, 'mtop_oc_sensitivity_convergence_time_compare');
close(fig);

end

function cmap = blue_white_red(n)

if nargin < 1
  n = 256;
end

bottom = [0 0.2 0.8];
middle = [1 1 1];
top = [0.8 0 0];
x = linspace(0, 1, n)';
cmap = zeros(n, 3);
lower = x <= 0.5;
cmap(lower, :) = bottom + (middle - bottom) .* (x(lower) / 0.5);
cmap(~lower, :) = middle + (top - middle) .* ((x(~lower) - 0.5) / 0.5);

end

function export_figure(fig, figuresDir, baseName)

pngPath = fullfile(figuresDir, [baseName '.png']);
pdfPath = fullfile(figuresDir, [baseName '.pdf']);
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

end
