function result = run_mtop_mma_heaviside(cfg, eta)
%RUN_MTOP_MMA_HEAVISIDE MTOP MMA Heaviside-projection study for assignment step 4.
%
% This runner combines MMA with density filtering and a Heaviside
% projection following Wang, Lazarov and Sigmund (2011), as implemented
% in the lecture Chapter 9 example ex4.m for ft = 3. The projection
% threshold eta and the multiresolution density mesh are used as in
% Nguyen et al. (2010). Continuation on the projection sharpness beta
% follows the same schedule as the lecture example: beta is doubled
% every iterbetamax = 50 iterations, starting from beta = 1 and stopping
% at betamax = 16.

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

if nargin < 2 || isempty(eta)
  eta = 0.5;
end

p = cfg.parameters;

problem = struct();
problem.id = sprintf('mtop_mma_heaviside_eta_%03i', round(100 * eta));
problem.feNelx = p.mtopNelx;
problem.feNely = p.mtopNely;
problem.densityPerElement = p.densityPerElement;
problem.densityNelx = p.mtopDensityNelx;
problem.densityNely = p.mtopDensityNely;
problem.volfrac = p.volfrac;
problem.penal = p.penal;
problem.rmin = p.classicalFilterRadius;
problem.tol = p.tol;
problem.maxIter = 800;
problem.E0 = 1;
problem.Emin = 1e-9;
problem.nu = 0.3;
problem.eta = eta;
problem.beta0 = p.heavisideBeta0;
problem.betaMax = p.heavisideBetaMax;
problem.iterBetaMax = 50;

fprintf('\nRunning assignment step 4: MTOP MMA Heaviside (eta = %.2f)\n', eta);
fprintf('  FE mesh:       %i x %i elements\n', problem.feNelx, problem.feNely);
fprintf('  Density mesh:  %i x %i cells\n', problem.densityNelx, problem.densityNely);
fprintf('  Cells per FE:  %i x %i\n', problem.densityPerElement(1), problem.densityPerElement(2));
fprintf('  Volume target: %.4f\n', problem.volfrac);
fprintf('  Penalization:  %.4f\n', problem.penal);
fprintf('  Filter radius: %.4f density cells\n', problem.rmin);
fprintf('  Heaviside eta: %.2f\n', problem.eta);
fprintf('  Beta range:    %.2f -> %.2f (doubled every %i iters)\n\n', ...
  problem.beta0, problem.betaMax, problem.iterBetaMax);

result = solve_mtop_mma_heaviside(problem);
result.problem = problem;

result.fineMeshCompliance = evaluate_fine_mbb_compliance(result.xPhys, problem);

resultPath = fullfile(cfg.paths.results, [problem.id '.mat']);
historyPath = fullfile(cfg.paths.results, [problem.id '_history.csv']);
summaryPath = fullfile(cfg.paths.results, [problem.id '_summary.txt']);

save(resultPath, 'result', '-v7.3');
write_history_csv(historyPath, result.history);
write_summary(summaryPath, result);
export_mma_design_figure(cfg.paths.figures, result);
export_mma_heaviside_convergence(cfg.paths.figures, result);

fprintf('\nSaved step 4 (MTOP MMA Heaviside, eta = %.2f) result data:\n', eta);
fprintf('  %s\n', resultPath);
fprintf('  %s\n', historyPath);
fprintf('  %s\n', summaryPath);
fprintf('Saved step 4 (MTOP MMA Heaviside, eta = %.2f) figures in:\n', eta);
fprintf('  %s\n\n', cfg.paths.figures);

end

function result = solve_mtop_mma_heaviside(problem)

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
maxIter = problem.maxIter;
E0 = problem.E0;
Emin = problem.Emin;
nu = problem.nu;
eta = problem.eta;
betaMax = problem.betaMax;
iterBetaMax = problem.iterBetaMax;

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
xmin = 0;
xmax = 1;

beta = problem.beta0;

history = struct();
history.iteration = zeros(maxIter, 1);
history.timeSeconds = zeros(maxIter, 1);
history.compliance = zeros(maxIter, 1);
history.volumeFraction = zeros(maxIter, 1);
history.change = zeros(maxIter, 1);
history.beta = zeros(maxIter, 1);
history.timeFEAssembly = zeros(maxIter, 1);
history.timeFESolve = zeros(maxIter, 1);
history.timeFilter = zeros(maxIter, 1);
history.timeOptimizer = zeros(maxIter, 1);

timer = tic;
change = 1;
iter = 0;
iterBeta = 1;
mmaparams = [];

while (change > tol || beta < betaMax) && iter < maxIter
  iter = iter + 1;

  tFilter = tic;
  xTilde = conv2(x, h, 'same') ./ Hs;
  xPhys = (tanh(beta * eta) + tanh(beta * (xTilde - eta))) ./ ...
    (tanh(beta * eta) + tanh(beta * (1 - eta)));
  filterTime0 = toc(tFilter);

  tAssembly = tic;
  sK = mtop_assemble_stiffness(xPhys, I_cells, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);
  K = sparse(iK, jK, sK);
  K = (K + K') / 2;
  assemblyTime = toc(tAssembly);

  tSolve = tic;
  U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
  solveTime = toc(tSolve);

  ceCells = mtop_strain_energy(U(edofMat), I_cells, ...
    densityPerX, densityPerY, feNelx, feNely);
  c = sum(sum((Emin + (E0 - Emin) * xPhys.^penal) .* ceCells));
  v = mean(xPhys(:));
  f0 = c / 100;
  f = v / volfrac - 1;

  tFilter = tic;
  dcdrho = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ceCells;
  dvdrho = ones(densityNely, densityNelx) / (densityNelx * densityNely);
  dxPhysdxTilde = beta * sech(beta * (xTilde - eta)).^2 ./ ...
    (tanh(beta * eta) + tanh(beta * (1 - eta)));
  dcdxTilde = dcdrho .* dxPhysdxTilde;
  dvdxTilde = dvdrho .* dxPhysdxTilde;
  dcdx = conv2(dcdxTilde ./ Hs, h, 'same');
  dvdx = conv2(dvdxTilde ./ Hs, h, 'same');
  df0dx = dcdx(:)' / 100;
  dfdx = dvdx(:)' / volfrac;
  filterTime1 = toc(tFilter);

  tOpt = tic;
  [xnew, ~, ~, ~, mmaparams, ~, change] = mma(x(:), xmin, xmax, f0, f, df0dx, dfdx, mmaparams, 'silent');
  xnew = reshape(xnew, densityNely, densityNelx);
  optimizerTime = toc(tOpt);

  elapsed = toc(timer);

  history.iteration(iter) = iter;
  history.timeSeconds(iter) = elapsed;
  history.compliance(iter) = c;
  history.volumeFraction(iter) = v;
  history.change(iter) = change;
  history.beta(iter) = beta;
  history.timeFEAssembly(iter) = assemblyTime;
  history.timeFESolve(iter) = solveTime;
  history.timeFilter(iter) = filterTime0 + filterTime1;
  history.timeOptimizer(iter) = optimizerTime;

  fprintf('[%s] iter: %4i | beta: %5.2f | change: %9.5f | c: %12.5f | v: %9.5f | time: %9.2f s\n', ...
    datestr(now, 'HH:MM:SS'), iter, beta, change, c, v, elapsed);

  if beta < betaMax && iterBeta >= iterBetaMax
    beta = min(2 * beta, betaMax);
    iterBeta = 0;
    change = 1;
    fprintf('  Heaviside beta increased to %.2f\n', beta);
  end

  if change > tol || beta < betaMax
    iterBeta = iterBeta + 1;
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
result.finalBeta = beta;
result.iterations = iter;
result.elapsedSeconds = history.timeSeconds(end);
result.converged = (change <= tol) && (beta >= betaMax);
result.filterKernel = h;
result.filterWeights = Hs;
result.timeBreakdown = struct();
result.timeBreakdown.feAssembly = sum(history.timeFEAssembly);
result.timeBreakdown.feSolve = sum(history.timeFESolve);
result.timeBreakdown.filter = sum(history.timeFilter);
result.timeBreakdown.optimizer = sum(history.timeOptimizer);
result.timeBreakdown.total = result.elapsedSeconds;

if ~result.converged
  warning('run_mtop_mma_heaviside:notConverged', ...
    'Maximum iteration count reached before the change tolerance was met at beta = betaMax.');
end

end


function sK = mtop_assemble_stiffness(xPhys, I_cells, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin)

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
  history.volumeFraction, history.change, history.beta, ...
  history.timeFEAssembly, history.timeFESolve, ...
  history.timeFilter, history.timeOptimizer, ...
  'VariableNames', {'iteration', 'time_seconds', 'compliance', 'volume_fraction', 'change', 'beta', ...
    'time_fe_assembly', 'time_fe_solve', 'time_filter', 'time_optimizer'});
writetable(T, historyPath);

end

function write_summary(summaryPath, result)

fid = fopen(summaryPath, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'MTOP MMA Heaviside (eta = %.4f)\n', result.problem.eta);
fprintf(fid, 'FE mesh: %i x %i elements\n', result.problem.feNelx, result.problem.feNely);
fprintf(fid, 'Density mesh: %i x %i cells\n', result.problem.densityNelx, result.problem.densityNely);
fprintf(fid, 'Density cells per FE: %i x %i\n', result.problem.densityPerElement(1), result.problem.densityPerElement(2));
fprintf(fid, 'Volume fraction target: %.6f\n', result.problem.volfrac);
fprintf(fid, 'Penalization power: %.6f\n', result.problem.penal);
fprintf(fid, 'Density filter radius: %.6f density cells\n', result.problem.rmin);
fprintf(fid, 'Heaviside eta: %.6f\n', result.problem.eta);
fprintf(fid, 'Heaviside beta range: %.4f -> %.4f\n', result.problem.beta0, result.problem.betaMax);
fprintf(fid, 'Iterations per beta level: %i\n', result.problem.iterBetaMax);
fprintf(fid, 'Change tolerance: %.6f\n', result.problem.tol);
fprintf(fid, 'Converged: %i\n', result.converged);
fprintf(fid, 'Iterations: %i\n', result.iterations);
fprintf(fid, 'Final beta: %.6f\n', result.finalBeta);
fprintf(fid, 'Elapsed time [s]: %.6f\n', result.elapsedSeconds);
fprintf(fid, 'Final MTOP compliance: %.12f\n', result.finalCompliance);
fprintf(fid, 'Final fine-mesh compliance: %.12f\n', result.fineMeshCompliance);
fprintf(fid, 'Final volume fraction: %.12f\n', result.finalVolumeFraction);
fprintf(fid, 'Final maximum design change: %.12f\n', result.finalChange);
fprintf(fid, 'Cumulative FE assembly time [s]: %.6f\n', result.timeBreakdown.feAssembly);
fprintf(fid, 'Cumulative FE solve time [s]: %.6f\n', result.timeBreakdown.feSolve);
fprintf(fid, 'Cumulative filter time [s]: %.6f\n', result.timeBreakdown.filter);
fprintf(fid, 'Cumulative optimizer time [s]: %.6f\n', result.timeBreakdown.optimizer);

end

function export_mma_design_figure(figuresDir, result)

if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(1 - result.xPhys);
colormap(fig, gray);
caxis([0 1]);
axis equal tight off;
export_figure(fig, figuresDir, [result.id '_design']);
close(fig);

end

function export_mma_heaviside_convergence(figuresDir, result)

history = result.history;

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
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
nexttile;
stairs(history.iteration, history.beta, 'LineWidth', 1.4);
grid on;
xlabel('Iteration number');
ylabel('Heaviside \beta');
title('Heaviside continuation schedule');
export_figure(fig, figuresDir, [result.id '_convergence_iter']);
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
export_figure(fig, figuresDir, [result.id '_convergence_time']);
close(fig);

end

function export_figure(fig, figuresDir, baseName)

pngPath = fullfile(figuresDir, [baseName '.png']);
pdfPath = fullfile(figuresDir, [baseName '.pdf']);
exportgraphics(fig, pngPath, 'Resolution', 300);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

end
