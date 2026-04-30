function result = run_classical_oc_density(cfg)
%RUN_CLASSICAL_OC_DENSITY Classical MBB OC baseline with density filtering for assignment step 3.
%
% This runner mirrors run_classical_oc_sensitivity but replaces the
% sensitivity filter by a density filter, in which the design variables x
% are convolved with the cone kernel H to obtain the physical densities
% rho on which the SIMP analysis is performed. The objective and constraint
% sensitivities with respect to the design variables follow from the chain
% rule applied to the linear filter operator.

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

p = cfg.parameters;

problem = struct();
problem.id = 'classical_oc_density';
problem.nelx = p.classicalNelx;
problem.nely = p.classicalNely;
problem.volfrac = p.volfrac;
problem.penal = p.penal;
problem.rmin = p.classicalFilterRadius;
problem.tol = p.tol;
problem.move = p.move;
problem.maxIter = 300;
problem.E0 = 1;
problem.Emin = 1e-9;
problem.nu = 0.3;

fprintf('\nRunning assignment step 3: classical OC density baseline\n');
fprintf('  Mesh:          %i x %i elements\n', problem.nelx, problem.nely);
fprintf('  Volume target: %.4f\n', problem.volfrac);
fprintf('  Penalization:  %.4f\n', problem.penal);
fprintf('  Filter radius: %.4f elements (density filter)\n\n', problem.rmin);

result = solve_classical_mbb_oc_density(problem);
result.problem = problem;

resultPath = fullfile(cfg.paths.results, [problem.id '.mat']);
historyPath = fullfile(cfg.paths.results, [problem.id '_history.csv']);
summaryPath = fullfile(cfg.paths.results, [problem.id '_summary.txt']);

save(resultPath, 'result', '-v7.3');
write_history_csv(historyPath, result.history);
write_summary(summaryPath, result);
export_classical_oc_figures(cfg.paths.figures, result);

fprintf('\nSaved step 3 (classical) result data:\n');
fprintf('  %s\n', resultPath);
fprintf('  %s\n', historyPath);
fprintf('  %s\n', summaryPath);
fprintf('Saved step 3 (classical) figures in:\n');
fprintf('  %s\n\n', cfg.paths.figures);

end

function result = solve_classical_mbb_oc_density(problem)

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

[dy, dx] = meshgrid(-ceil(rmin) + 1:ceil(rmin) - 1, -ceil(rmin) + 1:ceil(rmin) - 1);
h = max(0, rmin - sqrt(dx.^2 + dy.^2));
Hs = conv2(ones(nely, nelx), h, 'same');

x = repmat(volfrac, nely, nelx);
xPhys = conv2(x, h, 'same') ./ Hs;

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

  xPhys = conv2(x, h, 'same') ./ Hs;

  sK = reshape(KE(:) * (Emin + xPhys(:)'.^penal * (E0 - Emin)), 64 * nelx * nely, 1);
  K = sparse(iK, jK, sK);
  K = (K + K') / 2;
  U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

  ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
  c = sum(sum((Emin + xPhys.^penal * (E0 - Emin)) .* ce));
  v = sum(xPhys(:)) / (nelx * nely);

  dcdxPhys = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ce;
  dvdxPhys = ones(nely, nelx) / (nelx * nely);
  dcdx = conv2(dcdxPhys ./ Hs, h, 'same');
  dvdx = conv2(dvdxPhys ./ Hs, h, 'same');

  l1 = 0;
  l2 = 1e9;
  while (l2 - l1) / (l1 + l2) > 1e-3
    lmid = 0.5 * (l2 + l1);
    xCandidate = x .* sqrt(-dcdx ./ dvdx / lmid);
    xnew = max(0, max(x - move, min(1, min(x + move, xCandidate))));
    xnewPhys = conv2(xnew, h, 'same') ./ Hs;
    if sum(xnewPhys(:)) > volfrac * nelx * nely
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
  warning('run_classical_oc_density:notConverged', ...
    'Maximum iteration count reached before the change tolerance was met.');
end

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

fprintf(fid, 'Classical OC density baseline\n');
fprintf(fid, 'Mesh: %i x %i elements\n', result.problem.nelx, result.problem.nely);
fprintf(fid, 'Volume fraction target: %.6f\n', result.problem.volfrac);
fprintf(fid, 'Penalization power: %.6f\n', result.problem.penal);
fprintf(fid, 'Density filter radius: %.6f elements\n', result.problem.rmin);
fprintf(fid, 'Move limit: %.6f\n', result.problem.move);
fprintf(fid, 'Change tolerance: %.6f\n', result.problem.tol);
fprintf(fid, 'Converged: %i\n', result.converged);
fprintf(fid, 'Iterations: %i\n', result.iterations);
fprintf(fid, 'Elapsed time [s]: %.6f\n', result.elapsedSeconds);
fprintf(fid, 'Final compliance: %.12f\n', result.finalCompliance);
fprintf(fid, 'Final volume fraction: %.12f\n', result.finalVolumeFraction);
fprintf(fid, 'Final maximum design change: %.12f\n', result.finalChange);

end

function export_classical_oc_figures(figuresDir, result)

if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

history = result.history;
baseId = result.id;

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(1 - result.xPhys);
colormap(fig, gray);
caxis([0 1]);
axis equal tight off;
export_figure(fig, figuresDir, [baseId '_design']);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(history.iteration, history.compliance, 'LineWidth', 1.4);
grid on;
xlabel('Iteration number');
ylabel('Compliance');
title('Compliance convergence by iteration');
nexttile;
plot(history.iteration, history.volumeFraction, 'LineWidth', 1.4);
hold on;
yline(result.problem.volfrac, ':', 'Target', 'LineWidth', 1.0);
grid on;
xlabel('Iteration number');
ylabel('Volume fraction');
title('Volume convergence by iteration');
export_figure(fig, figuresDir, [baseId '_convergence_iter']);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(history.timeSeconds, history.compliance, 'LineWidth', 1.4);
grid on;
xlabel('Computation time [s]');
ylabel('Compliance');
title('Compliance convergence by time');
nexttile;
plot(history.timeSeconds, history.volumeFraction, 'LineWidth', 1.4);
hold on;
yline(result.problem.volfrac, ':', 'Target', 'LineWidth', 1.0);
grid on;
xlabel('Computation time [s]');
ylabel('Volume fraction');
title('Volume convergence by time');
export_figure(fig, figuresDir, [baseId '_convergence_time']);
close(fig);

end

function export_figure(fig, figuresDir, baseName)

pngPath = fullfile(figuresDir, [baseName '.png']);
pdfPath = fullfile(figuresDir, [baseName '.pdf']);
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

end
