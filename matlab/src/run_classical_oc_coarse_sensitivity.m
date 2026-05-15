function result = run_classical_oc_coarse_sensitivity(cfg)

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

p = cfg.parameters;

problem = struct();
problem.id = 'classical_oc_coarse_sensitivity';
problem.nelx = p.mtopNelx;
problem.nely = p.mtopNely;
problem.fineNelx = p.classicalNelx;
problem.fineNely = p.classicalNely;
problem.volfrac = p.volfrac;
problem.penal = p.penal;
problem.rmin = p.classicalFilterRadius / p.densityPerElement(1);
problem.tol = p.tol;
problem.move = p.move;
problem.maxIter = 300;
problem.E0 = 1;
problem.Emin = 1e-9;
problem.nu = 0.3;

fprintf('\nRunning coarse-coarse OC sensitivity baseline\n');
fprintf('  Coarse mesh:   %i x %i elements\n', problem.nelx, problem.nely);
fprintf('  Fine check:    %i x %i elements\n', problem.fineNelx, problem.fineNely);
fprintf('  Volume target: %.4f\n', problem.volfrac);
fprintf('  Penalization:  %.4f\n', problem.penal);
fprintf('  Filter radius: %.4f coarse elements\n\n', problem.rmin);

result = solve_classical_mbb_oc(problem);
result.problem = problem;
result.xFine = upscale_density(result.xPhys, problem);
result.fineMeshCompliance = evaluate_fine_mbb_compliance(result.xFine, problem);
result.grayness = mean(4 * result.xPhys(:) .* (1 - result.xPhys(:)));

resultPath = fullfile(cfg.paths.results, [problem.id '.mat']);
historyPath = fullfile(cfg.paths.results, [problem.id '_history.csv']);
summaryPath = fullfile(cfg.paths.results, [problem.id '_summary.txt']);

save(resultPath, 'result', '-v7.3');
write_history_csv(historyPath, result.history);
write_summary(summaryPath, result);
export_coarse_figures(cfg.paths.figures, result);

fprintf('\nSaved coarse-coarse baseline data:\n');
fprintf('  %s\n', resultPath);
fprintf('  %s\n', historyPath);
fprintf('  %s\n', summaryPath);
fprintf('Saved coarse-coarse figures in:\n');
fprintf('  %s\n\n', cfg.paths.figures);

end

function result = solve_classical_mbb_oc(problem)

nelx = problem.nelx;
nely = problem.nely;
volfrac = problem.volfrac;
penal = problem.penal;
rmin = problem.rmin;
tol = problem.tol;
move = problem.move;
maxIter = problem.maxIter;
E0 = problem.E0;
Emin = problem.Emin;
nu = problem.nu;

KE = element_stiffness_matrix(nu);
[edofMat, iK, jK, F, freedofs] = mbb_topology(nelx, nely);
U = zeros(2 * (nely + 1) * (nelx + 1), 1);

[dy, dx] = meshgrid(-ceil(rmin) + 1:ceil(rmin) - 1, -ceil(rmin) + 1:ceil(rmin) - 1);
h = max(0, rmin - sqrt(dx.^2 + dy.^2));
Hs = conv2(ones(nely, nelx), h, 'same');

x = repmat(volfrac, nely, nelx);
xPhys = x;

history = struct();
history.iteration = zeros(maxIter, 1);
history.timeSeconds = zeros(maxIter, 1);
history.compliance = zeros(maxIter, 1);
history.volumeFraction = zeros(maxIter, 1);
history.change = zeros(maxIter, 1);
history.timeFEAssembly = zeros(maxIter, 1);
history.timeFESolve = zeros(maxIter, 1);
history.timeFilter = zeros(maxIter, 1);
history.timeOptimizer = zeros(maxIter, 1);

timer = tic;
change = 1;
iter = 0;

while change > tol && iter < maxIter
  iter = iter + 1;
  xPhys = x;

  tAssembly = tic;
  sK = reshape(KE(:) * (Emin + xPhys(:)'.^penal * (E0 - Emin)), 64 * nelx * nely, 1);
  K = sparse(iK, jK, sK);
  K = (K + K') / 2;
  assemblyTime = toc(tAssembly);

  tSolve = tic;
  U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
  solveTime = toc(tSolve);

  ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
  c = sum(sum((Emin + xPhys.^penal * (E0 - Emin)) .* ce));
  v = mean(xPhys(:));

  tFilter = tic;
  dcdx = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ce;
  dvdx = ones(nely, nelx) / (nelx * nely);
  dcdx = conv2(dcdx .* xPhys, h, 'same') ./ Hs ./ max(1e-3, xPhys);
  filterTime = toc(tFilter);

  tOpt = tic;
  l1 = 0;
  l2 = 1e9;
  while (l2 - l1) / (l1 + l2) > 1e-3
    lmid = 0.5 * (l2 + l1);
    xCandidate = x .* sqrt(-dcdx ./ dvdx / lmid);
    xnew = max(0, max(x - move, min(1, min(x + move, xCandidate))));
    if sum(xnew(:)) > volfrac * nelx * nely
      l1 = lmid;
    else
      l2 = lmid;
    end
  end
  optimizerTime = toc(tOpt);

  change = max(abs(xnew(:) - x(:)));
  elapsed = toc(timer);

  history.iteration(iter) = iter;
  history.timeSeconds(iter) = elapsed;
  history.compliance(iter) = c;
  history.volumeFraction(iter) = v;
  history.change(iter) = change;
  history.timeFEAssembly(iter) = assemblyTime;
  history.timeFESolve(iter) = solveTime;
  history.timeFilter(iter) = filterTime;
  history.timeOptimizer(iter) = optimizerTime;

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
result.timeBreakdown = struct();
result.timeBreakdown.feAssembly = sum(history.timeFEAssembly);
result.timeBreakdown.feSolve = sum(history.timeFESolve);
result.timeBreakdown.filter = sum(history.timeFilter);
result.timeBreakdown.optimizer = sum(history.timeOptimizer);
result.timeBreakdown.total = result.elapsedSeconds;

end

function xFine = upscale_density(xCoarse, problem)

scaleX = problem.fineNelx / problem.nelx;
scaleY = problem.fineNely / problem.nely;
if abs(scaleX - round(scaleX)) > eps || abs(scaleY - round(scaleY)) > eps
  error('run_classical_oc_coarse_sensitivity:nonIntegerScale', ...
    'Fine mesh must be an integer refinement of the coarse mesh.');
end
xFine = kron(xCoarse, ones(round(scaleY), round(scaleX)));

end

function compliance = evaluate_fine_mbb_compliance(xPhys, problem)

nelx = size(xPhys, 2);
nely = size(xPhys, 1);
KE = element_stiffness_matrix(problem.nu);
[edofMat, iK, jK, F, freedofs] = mbb_topology(nelx, nely);

sK = reshape(KE(:) * (problem.Emin + xPhys(:)'.^problem.penal * ...
  (problem.E0 - problem.Emin)), 64 * nelx * nely, 1);
K = sparse(iK, jK, sK);
K = (K + K') / 2;
U = zeros(2 * (nely + 1) * (nelx + 1), 1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
compliance = sum(sum((problem.Emin + xPhys.^problem.penal * ...
  (problem.E0 - problem.Emin)) .* ce));

end

function [edofMat, iK, jK, F, freedofs] = mbb_topology(nelx, nely)

nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(2 * nodenrs(1:end-1, 1:end-1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1);
F = sparse(2, 1, -1, 2 * (nely + 1) * (nelx + 1), 1);
fixeddofs = union(1:2:2 * (nely + 1), 2 * (nelx + 1) * (nely + 1));
alldofs = 1:2 * (nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

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
  history.timeFEAssembly, history.timeFESolve, ...
  history.timeFilter, history.timeOptimizer, ...
  'VariableNames', {'iteration', 'time_seconds', 'compliance', 'volume_fraction', ...
  'change', 'time_fe_assembly', 'time_fe_solve', 'time_filter', 'time_optimizer'});
writetable(T, historyPath);

end

function write_summary(summaryPath, result)

fid = fopen(summaryPath, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Classical OC coarse-coarse sensitivity baseline\n');
fprintf(fid, 'Coarse mesh: %i x %i elements\n', result.problem.nelx, result.problem.nely);
fprintf(fid, 'Fine reanalysis mesh: %i x %i elements\n', result.problem.fineNelx, result.problem.fineNely);
fprintf(fid, 'Volume fraction target: %.6f\n', result.problem.volfrac);
fprintf(fid, 'Penalization power: %.6f\n', result.problem.penal);
fprintf(fid, 'Sensitivity filter radius: %.6f coarse elements\n', result.problem.rmin);
fprintf(fid, 'Move limit: %.6f\n', result.problem.move);
fprintf(fid, 'Change tolerance: %.6f\n', result.problem.tol);
fprintf(fid, 'Converged: %i\n', result.converged);
fprintf(fid, 'Iterations: %i\n', result.iterations);
fprintf(fid, 'Elapsed time [s]: %.6f\n', result.elapsedSeconds);
fprintf(fid, 'Native coarse compliance: %.12f\n', result.finalCompliance);
fprintf(fid, 'Fine-mesh compliance of upsampled coarse design: %.12f\n', result.fineMeshCompliance);
fprintf(fid, 'Final volume fraction: %.12f\n', result.finalVolumeFraction);
fprintf(fid, 'Final maximum design change: %.12f\n', result.finalChange);
fprintf(fid, 'Grayness index mean(4*rho*(1-rho)): %.12f\n', result.grayness);
fprintf(fid, 'Cumulative FE assembly time [s]: %.6f\n', result.timeBreakdown.feAssembly);
fprintf(fid, 'Cumulative FE solve time [s]: %.6f\n', result.timeBreakdown.feSolve);
fprintf(fid, 'Cumulative filter time [s]: %.6f\n', result.timeBreakdown.filter);
fprintf(fid, 'Cumulative optimizer time [s]: %.6f\n', result.timeBreakdown.optimizer);

end

function export_coarse_figures(figuresDir, result)

if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(1 - result.xFine);
colormap(fig, gray);
caxis([0 1]);
axis equal tight off;
export_figure(fig, figuresDir, 'classical_oc_coarse_sensitivity_design');
close(fig);

history = result.history;
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(history.iteration, history.compliance, 'LineWidth', 1.4);
grid on;
xlabel('Iteration number');
ylabel('Compliance');
title('Coarse-coarse compliance convergence');
nexttile;
plot(history.iteration, history.volumeFraction, 'LineWidth', 1.4);
hold on;
yline(result.problem.volfrac, ':', 'Target', 'LineWidth', 1.0);
grid on;
xlabel('Iteration number');
ylabel('Volume fraction');
title('Coarse-coarse volume convergence');
export_figure(fig, figuresDir, 'classical_oc_coarse_sensitivity_convergence_iter');
close(fig);

end

function export_figure(fig, figuresDir, baseName)

exportgraphics(fig, fullfile(figuresDir, [baseName '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(figuresDir, [baseName '.pdf']), 'ContentType', 'vector');

end
