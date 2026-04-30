function report = verify_sensitivities(cfg, options)
%VERIFY_SENSITIVITIES Finite-difference verification of compliance and volume sensitivities.
%
% This routine compares the analytical sensitivities used by the assignment
% codes with central finite-difference approximations. The verification is
% performed on a small instance of the MBB problem because the FD evaluation
% requires repeated FE solves. The same SIMP, BC and integration conventions
% are used as in run_classical_oc_sensitivity and run_mtop_oc_sensitivity,
% so the test exercises the actual derivatives that drive the optimization.
%
% A successful test produces a relative error below ~1e-5 for compliance
% and below ~1e-12 for the (linear) volume sensitivity.

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

if nargin < 2 || isempty(options)
  options = struct();
end
if ~isfield(options, 'feNelx'), options.feNelx = 12; end
if ~isfield(options, 'feNely'), options.feNely = 4; end
if ~isfield(options, 'densityPerX'), options.densityPerX = cfg.parameters.densityPerElement(1); end
if ~isfield(options, 'densityPerY'), options.densityPerY = cfg.parameters.densityPerElement(2); end
if ~isfield(options, 'numSamples'), options.numSamples = 20; end
if ~isfield(options, 'fdStep'), options.fdStep = 1e-6; end
if ~isfield(options, 'penal'), options.penal = cfg.parameters.penal; end
if ~isfield(options, 'volfrac'), options.volfrac = cfg.parameters.volfrac; end
if ~isfield(options, 'E0'), options.E0 = 1; end
if ~isfield(options, 'Emin'), options.Emin = 1e-9; end
if ~isfield(options, 'nu'), options.nu = 0.3; end
if ~isfield(options, 'seed'), options.seed = 0; end

reportClassical = check_classical(options);
reportMtop = check_mtop(options);
reportHeaviside = check_heaviside_chain(options);

report = struct();
report.classical = reportClassical;
report.mtop = reportMtop;
report.heaviside = reportHeaviside;

print_report(report);

if isfield(cfg, 'paths') && isfield(cfg.paths, 'results')
  outputPath = fullfile(cfg.paths.results, 'sensitivity_verification.mat');
  save(outputPath, 'report', '-v7.3');
  fprintf('Saved sensitivity verification report:\n  %s\n', outputPath);
end

end


function out = check_classical(options)

rng(options.seed);

nelx = options.feNelx;
nely = options.feNely;
penal = options.penal;
E0 = options.E0;
Emin = options.Emin;
nu = options.nu;

KE = element_stiffness_matrix(nu);
[edofMat, iK, jK, F, freedofs] = mbb_topology(nelx, nely);

xPhys = options.volfrac + 0.05 * (rand(nely, nelx) - 0.5);
xPhys = min(max(xPhys, 0.05), 0.95);

[c0, dcdrho, dvdrho] = classical_objective(xPhys, KE, edofMat, iK, jK, F, freedofs, ...
  penal, E0, Emin);

numSamples = min(options.numSamples, numel(xPhys));
sampleIdx = randperm(numel(xPhys), numSamples);

dcFD = zeros(numSamples, 1);
dvFD = zeros(numSamples, 1);
for k = 1:numSamples
  i = sampleIdx(k);
  xPlus = xPhys; xPlus(i) = xPlus(i) + options.fdStep;
  xMinus = xPhys; xMinus(i) = xMinus(i) - options.fdStep;
  cPlus  = classical_objective(xPlus,  KE, edofMat, iK, jK, F, freedofs, penal, E0, Emin);
  cMinus = classical_objective(xMinus, KE, edofMat, iK, jK, F, freedofs, penal, E0, Emin);
  dcFD(k) = (cPlus - cMinus) / (2 * options.fdStep);
  vPlus  = mean(xPlus(:));
  vMinus = mean(xMinus(:));
  dvFD(k) = (vPlus - vMinus) / (2 * options.fdStep);
end

dcAna = dcdrho(sampleIdx);
dvAna = dvdrho(sampleIdx);

out = error_metrics('classical', c0, sampleIdx, dcAna, dcFD, dvAna, dvFD);

end


function out = check_mtop(options)

rng(options.seed);

feNelx = options.feNelx;
feNely = options.feNely;
densityPerX = options.densityPerX;
densityPerY = options.densityPerY;
penal = options.penal;
E0 = options.E0;
Emin = options.Emin;
nu = options.nu;

I_cells = mtop_subcell_stiffness(densityPerX, densityPerY, nu);
[edofMat, iK, jK, F, freedofs] = mbb_topology(feNelx, feNely);

densityNelx = feNelx * densityPerX;
densityNely = feNely * densityPerY;
xPhys = options.volfrac + 0.05 * (rand(densityNely, densityNelx) - 0.5);
xPhys = min(max(xPhys, 0.05), 0.95);

[c0, dcdrho, dvdrho] = mtop_objective(xPhys, I_cells, edofMat, iK, jK, F, freedofs, ...
  densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);

numSamples = min(options.numSamples, numel(xPhys));
sampleIdx = randperm(numel(xPhys), numSamples);

dcFD = zeros(numSamples, 1);
dvFD = zeros(numSamples, 1);
for k = 1:numSamples
  i = sampleIdx(k);
  xPlus = xPhys; xPlus(i) = xPlus(i) + options.fdStep;
  xMinus = xPhys; xMinus(i) = xMinus(i) - options.fdStep;
  cPlus  = mtop_objective(xPlus,  I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);
  cMinus = mtop_objective(xMinus, I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);
  dcFD(k) = (cPlus - cMinus) / (2 * options.fdStep);
  vPlus  = mean(xPlus(:));
  vMinus = mean(xMinus(:));
  dvFD(k) = (vPlus - vMinus) / (2 * options.fdStep);
end

dcAna = dcdrho(sampleIdx);
dvAna = dvdrho(sampleIdx);

out = error_metrics('mtop', c0, sampleIdx, dcAna, dcFD, dvAna, dvFD);

end


function out = check_heaviside_chain(options)
%CHECK_HEAVISIDE_CHAIN Verify the full filter + Heaviside + MTOP chain rule.
%
% Builds a small density grid, applies the cone filter, then the Heaviside
% projection with a representative beta and eta, then evaluates compliance
% via the multiresolution scheme. The analytical sensitivity dC/dx is
% compared with central finite differences applied directly to the design
% variables x. A passing test verifies the entire chain that drives the
% MTOP MMA Heaviside runner of step 4.

rng(options.seed);

feNelx = options.feNelx;
feNely = options.feNely;
densityPerX = options.densityPerX;
densityPerY = options.densityPerY;
penal = options.penal;
E0 = options.E0;
Emin = options.Emin;
nu = options.nu;

beta = 8;
eta = 0.5;
rmin = 6;

I_cells = mtop_subcell_stiffness(densityPerX, densityPerY, nu);
[edofMat, iK, jK, F, freedofs] = mbb_topology(feNelx, feNely);

densityNelx = feNelx * densityPerX;
densityNely = feNely * densityPerY;
[dy, dx] = meshgrid(-ceil(rmin) + 1:ceil(rmin) - 1, -ceil(rmin) + 1:ceil(rmin) - 1);
h = max(0, rmin - sqrt(dx.^2 + dy.^2));
Hs = conv2(ones(densityNely, densityNelx), h, 'same');

x0 = options.volfrac + 0.10 * (rand(densityNely, densityNelx) - 0.5);
x0 = min(max(x0, 0.05), 0.95);

[c0, dcdx] = heaviside_objective(x0, h, Hs, I_cells, edofMat, iK, jK, F, freedofs, ...
  densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin, beta, eta);

numSamples = min(options.numSamples, numel(x0));
sampleIdx = randperm(numel(x0), numSamples);

dcFD = zeros(numSamples, 1);
for k = 1:numSamples
  i = sampleIdx(k);
  xPlus = x0; xPlus(i) = xPlus(i) + options.fdStep;
  xMinus = x0; xMinus(i) = xMinus(i) - options.fdStep;
  cPlus  = heaviside_objective(xPlus,  h, Hs, I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin, beta, eta);
  cMinus = heaviside_objective(xMinus, h, Hs, I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin, beta, eta);
  dcFD(k) = (cPlus - cMinus) / (2 * options.fdStep);
end

dcAna = dcdx(sampleIdx);

out = struct();
out.label = 'heaviside';
out.compliance = c0;
out.beta = beta;
out.eta = eta;
out.filterRadius = rmin;
out.sampleIdx = sampleIdx(:);
out.dcAnalytical = dcAna(:);
out.dcFiniteDiff = dcFD;
out.dcMaxAbsError = max(abs(dcAna(:) - dcFD));
out.dcMaxRelError = max(abs(dcAna(:) - dcFD)) / max(max(abs(dcFD)), eps);

end


function [c, dcdx] = heaviside_objective(x, h, Hs, I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin, beta, eta)

xTilde = conv2(x, h, 'same') ./ Hs;
xPhys = (tanh(beta * eta) + tanh(beta * (xTilde - eta))) ./ ...
  (tanh(beta * eta) + tanh(beta * (1 - eta)));

sK = mtop_assemble_local(xPhys, I_cells, densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);
K = sparse(iK, jK, sK);
K = (K + K') / 2;
U = zeros(2 * (feNely + 1) * (feNelx + 1), 1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

ceCells = mtop_strain_energy_local(U(edofMat), I_cells, densityPerX, densityPerY, feNelx, feNely);
c = sum(sum((Emin + (E0 - Emin) * xPhys.^penal) .* ceCells));

if nargout > 1
  densityNelx = feNelx * densityPerX;
  densityNely = feNely * densityPerY;
  dcdrho = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ceCells;
  dxPhysdxTilde = beta * sech(beta * (xTilde - eta)).^2 ./ ...
    (tanh(beta * eta) + tanh(beta * (1 - eta)));
  dcdxTilde = dcdrho .* dxPhysdxTilde;
  dcdx = conv2(dcdxTilde ./ Hs, h, 'same');
  dcdx = reshape(dcdx, densityNely, densityNelx);
end

end


function sK = mtop_assemble_local(xPhys, I_cells, densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin)

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


function ceCells = mtop_strain_energy_local(ueK, I_cells, densityPerX, densityPerY, feNelx, feNely)

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


function [c, dcdrho, dvdrho] = classical_objective(xPhys, KE, edofMat, iK, jK, F, freedofs, ...
    penal, E0, Emin)

[nely, nelx] = size(xPhys);
sK = reshape(KE(:) * (Emin + xPhys(:)'.^penal * (E0 - Emin)), 64 * nelx * nely, 1);
K = sparse(iK, jK, sK);
K = (K + K') / 2;
U = zeros(2 * (nely + 1) * (nelx + 1), 1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
c = sum(sum((Emin + xPhys.^penal * (E0 - Emin)) .* ce));

if nargout > 1
  dcdrho = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ce;
  dvdrho = ones(nely, nelx) / (nelx * nely);
end

end


function [c, dcdrho, dvdrho] = mtop_objective(xPhys, I_cells, edofMat, iK, jK, F, freedofs, ...
    densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin)

sK = mtop_assemble(xPhys, I_cells, densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin);
K = sparse(iK, jK, sK);
K = (K + K') / 2;
U = zeros(2 * (feNely + 1) * (feNelx + 1), 1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

ceCells = strain_energy(U(edofMat), I_cells, densityPerX, densityPerY, feNelx, feNely);
c = sum(sum((Emin + (E0 - Emin) * xPhys.^penal) .* ceCells));

if nargout > 1
  dcdrho = -penal * (E0 - Emin) * xPhys.^(penal - 1) .* ceCells;
  densityNelx = densityPerX * feNelx;
  densityNely = densityPerY * feNely;
  dvdrho = ones(densityNely, densityNelx) / (densityNelx * densityNely);
end

end


function sK = mtop_assemble(xPhys, I_cells, densityPerX, densityPerY, feNelx, feNely, penal, E0, Emin)

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


function ceCells = strain_energy(ueK, I_cells, densityPerX, densityPerY, feNelx, feNely)

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


function KE = element_stiffness_matrix(nu)

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1 / (1 - nu^2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);

end


function out = error_metrics(label, c0, idx, dcAna, dcFD, dvAna, dvFD)

dcAna = dcAna(:);
dcFD = dcFD(:);
dvAna = dvAna(:);
dvFD = dvFD(:);

cScale = max(abs(dcFD));
out = struct();
out.label = label;
out.compliance = c0;
out.sampleIdx = idx(:);
out.dcAnalytical = dcAna;
out.dcFiniteDiff = dcFD;
out.dcMaxAbsError = max(abs(dcAna - dcFD));
out.dcMaxRelError = max(abs(dcAna - dcFD)) / max(cScale, eps);
out.dvAnalytical = dvAna;
out.dvFiniteDiff = dvFD;
out.dvMaxAbsError = max(abs(dvAna - dvFD));

end


function print_report(report)

fprintf('\nSensitivity verification\n');
fprintf('  Classical OC formulation\n');
fprintf('    Compliance:                       %.6f\n', report.classical.compliance);
fprintf('    Max abs FD error in dC/drho:      %.3e\n', report.classical.dcMaxAbsError);
fprintf('    Max rel FD error in dC/drho:      %.3e\n', report.classical.dcMaxRelError);
fprintf('    Max abs FD error in dV/drho:      %.3e\n', report.classical.dvMaxAbsError);
fprintf('  MTOP formulation (per-cell sensitivity)\n');
fprintf('    Compliance:                       %.6f\n', report.mtop.compliance);
fprintf('    Max abs FD error in dC/drho:      %.3e\n', report.mtop.dcMaxAbsError);
fprintf('    Max rel FD error in dC/drho:      %.3e\n', report.mtop.dcMaxRelError);
fprintf('    Max abs FD error in dV/drho:      %.3e\n', report.mtop.dvMaxAbsError);
fprintf('  Heaviside chain (filter + projection + MTOP)\n');
fprintf('    beta = %.2f, eta = %.2f\n', report.heaviside.beta, report.heaviside.eta);
fprintf('    Compliance:                       %.6f\n', report.heaviside.compliance);
fprintf('    Max abs FD error in dC/dx:        %.3e\n', report.heaviside.dcMaxAbsError);
fprintf('    Max rel FD error in dC/dx:        %.3e\n', report.heaviside.dcMaxRelError);

end
