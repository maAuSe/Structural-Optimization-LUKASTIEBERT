function export_density_sensitivity_figures(cfg)
%EXPORT_DENSITY_SENSITIVITY_FIGURES Visualize density-filter sensitivities.

if nargin < 1 || isempty(cfg)
  paths = setup_project();
  cfg = assignment_config(paths);
end

classicalData = load(fullfile(cfg.paths.results, 'classical_oc_density.mat'), 'result');
mtopData = load(fullfile(cfg.paths.results, 'mtop_oc_density.mat'), 'result');

[classicalRaw, classicalFiltered] = classical_density_sensitivities(classicalData.result);
[mtopRaw, mtopFiltered] = mtop_density_sensitivities(mtopData.result);

panels = {
  classicalRaw,      classicalData.result.xPhys, 'Classical raw sensitivity'
  mtopRaw,           mtopData.result.xPhys,      'MTOP raw sensitivity'
  classicalFiltered, classicalData.result.xPhys, 'Classical after density-filter chain rule'
  mtopFiltered,      mtopData.result.xPhys,      'MTOP after density-filter chain rule'
};

for k = 1:size(panels, 1)
  field = max(panels{k, 1}, 0);
  scale = percentile_value(field(:), 99);
  if ~(scale > 0)
    scale = max(field(:));
  end
  panels{k, 1} = min(field / max(scale, eps), 1);
end

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
  'Position', [2 2 18.4 6.4]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'tight');

for k = 1:size(panels, 1)
  ax = nexttile(layout, k);
  imagesc(ax, panels{k, 1});
  hold(ax, 'on');
  contour(ax, panels{k, 2}, [0.5 0.5], '-', 'LineWidth', 0.5, ...
    'Color', [0.85 0.1 0.1]);
  hold(ax, 'off');
  axis(ax, 'equal', 'tight', 'off');
  colormap(ax, parula(256));
  caxis(ax, [0 1]);
  title(ax, panels{k, 3}, 'FontName', 'Arial', 'FontSize', 8.0, ...
    'FontWeight', 'bold');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontName = 'Arial';
cb.FontSize = 7.2;
cb.Label.String = 'normalized sensitivity (clipped at 99th percentile)';

exportgraphics(fig, fullfile(cfg.paths.figures, 'density_filter_sensitivity_fields.png'), ...
  'Resolution', 300);
exportgraphics(fig, fullfile(cfg.paths.figures, 'density_filter_sensitivity_fields.pdf'), ...
  'ContentType', 'vector');
close(fig);

end


function [rawMagnitude, filteredMagnitude] = classical_density_sensitivities(result)

xPhys = result.xPhys;
[nely, nelx] = size(xPhys);
problem = result.problem;

KE = element_stiffness_matrix(problem.nu);
[edofMat, ~, ~, ~, ~] = mbb_topology(nelx, nely);

ce = reshape(sum((result.U(edofMat) * KE) .* result.U(edofMat), 2), nely, nelx);
dcdrho = -problem.penal * (problem.E0 - problem.Emin) * xPhys.^(problem.penal - 1) .* ce;
dcdx = conv2(dcdrho ./ result.filterWeights, result.filterKernel, 'same');

rawMagnitude = -dcdrho;
filteredMagnitude = -dcdx;

end


function [rawMagnitude, filteredMagnitude] = mtop_density_sensitivities(result)

problem = result.problem;
feNelx = problem.feNelx;
feNely = problem.feNely;
densityPerX = problem.densityPerElement(1);
densityPerY = problem.densityPerElement(2);

I_cells = mtop_subcell_stiffness(densityPerX, densityPerY, problem.nu);
[edofMat, ~, ~, ~, ~] = mbb_topology(feNelx, feNely);

ceCells = mtop_strain_energy(result.U(edofMat), I_cells, densityPerX, ...
  densityPerY, feNelx, feNely);
dcdrho = -problem.penal * (problem.E0 - problem.Emin) * ...
  result.xPhys.^(problem.penal - 1) .* ceCells;
dcdx = conv2(dcdrho ./ result.filterWeights, result.filterKernel, 'same');

rawMagnitude = -dcdrho;
filteredMagnitude = -dcdx;

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


function ceCells = mtop_strain_energy(ueK, I_cells, densityPerX, densityPerY, feNelx, feNely)

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


function value = percentile_value(values, percent)

values = sort(values(:));
if isempty(values)
  value = NaN;
  return;
end

position = 1 + (numel(values) - 1) * percent / 100;
lower = floor(position);
upper = ceil(position);
if lower == upper
  value = values(lower);
else
  value = values(lower) + (position - lower) * (values(upper) - values(lower));
end

end
