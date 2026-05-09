function export_cone_kernel_figure()
%EXPORT_CONE_KERNEL_FIGURE Create the report figure for the cone filter.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
figuresDir = fullfile(repoRoot, 'figures');
if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

rmin = 24;
offsets = -rmin:rmin;
[dx, dy] = meshgrid(offsets, offsets);
distance = sqrt(dx.^2 + dy.^2);
H = max(0, rmin - distance);
Hn = H / rmin;

theta = linspace(0, 2 * pi, 300);
profileDistance = linspace(0, 1.15 * rmin, 500);
profileWeight = max(0, rmin - profileDistance) / rmin;

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
  'Position', [2 2 17.5 10.5]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axSurface = nexttile(layout, 1, [2 1]);
surf(axSurface, dx, dy, Hn, Hn, 'EdgeColor', [0.72 0.72 0.72], ...
  'EdgeAlpha', 0.30, 'FaceAlpha', 0.98);
hold(axSurface, 'on');
plot3(axSurface, rmin * cos(theta), rmin * sin(theta), ...
  zeros(size(theta)), 'Color', [0.12 0.12 0.12], 'LineWidth', 1.0);
plot3(axSurface, 0, 0, 1, 'o', 'MarkerSize', 5.5, ...
  'MarkerFaceColor', [0.89 0.10 0.11], 'MarkerEdgeColor', 'w', ...
  'LineWidth', 0.8);
view(axSurface, -37, 28);
axis(axSurface, 'tight');
box(axSurface, 'on');
grid(axSurface, 'on');
xlabel(axSurface, '\Delta x [cells]');
ylabel(axSurface, '\Delta y [cells]');
zlabel(axSurface, 'H_{ij}/r_{min}');
title(axSurface, 'Cone weights around cell i');
zlim(axSurface, [0 1]);
colormap(axSurface, parula(256));

axStencil = nexttile(layout, 2);
imagesc(axStencil, offsets, offsets, Hn);
set(axStencil, 'YDir', 'normal');
hold(axStencil, 'on');
plot(axStencil, rmin * cos(theta), rmin * sin(theta), 'w-', 'LineWidth', 1.1);
plot(axStencil, 0, 0, 'o', 'MarkerSize', 5.5, ...
  'MarkerFaceColor', [0.89 0.10 0.11], 'MarkerEdgeColor', 'w', ...
  'LineWidth', 0.8);
axis(axStencil, 'image');
xlim(axStencil, [-rmin rmin]);
ylim(axStencil, [-rmin rmin]);
xlabel(axStencil, '\Delta x [cells]');
ylabel(axStencil, '\Delta y [cells]');
title(axStencil, 'Filter footprint');
cb = colorbar(axStencil);
cb.Label.String = 'normalized weight';

axProfile = nexttile(layout, 4);
area(axProfile, profileDistance, profileWeight, ...
  'FaceColor', [0.30 0.69 0.76], 'FaceAlpha', 0.28, ...
  'EdgeColor', 'none');
hold(axProfile, 'on');
plot(axProfile, profileDistance, profileWeight, ...
  'Color', [0.07 0.37 0.42], 'LineWidth', 2.2);
xline(axProfile, rmin, '--', 'Color', [0.18 0.18 0.18], ...
  'LineWidth', 1.0);
text(axProfile, 14.7, 0.13, 'r_{min} = 24 cells', ...
  'FontSize', 8.7, 'Color', [0.12 0.12 0.12], ...
  'BackgroundColor', 'w', 'Margin', 2);
grid(axProfile, 'on');
box(axProfile, 'on');
xlabel(axProfile, 'distance from cell i [cells]');
ylabel(axProfile, 'H_{ij}/r_{min}');
title(axProfile, 'Linear decay to zero');
xlim(axProfile, [0 1.15 * rmin]);
ylim(axProfile, [0 1.04]);

axesList = [axSurface axStencil axProfile];
set(axesList, 'FontName', 'Arial', 'FontSize', 8.8, 'LineWidth', 0.75);
set([axSurface axProfile], 'GridAlpha', 0.18, 'MinorGridAlpha', 0.10);

exportgraphics(fig, fullfile(figuresDir, 'cone_kernel_radius_24.png'), ...
  'Resolution', 300);
exportgraphics(fig, fullfile(figuresDir, 'cone_kernel_radius_24.pdf'), ...
  'ContentType', 'vector');
close(fig);
end
