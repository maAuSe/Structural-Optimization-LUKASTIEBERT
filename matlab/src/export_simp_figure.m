function export_simp_figure()

repoRoot = fileparts(fileparts(mfilename('fullpath')));
figuresDir = fullfile(repoRoot, 'figures');
if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

rho = linspace(0, 1, 600);
E0 = 1;
Emin = 1e-9;
penalties = [1 2 3 4];
colors = [
  0.18 0.49 0.72
  0.47 0.67 0.19
  0.85 0.33 0.10
  0.49 0.18 0.56
];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
  'Position', [2 2 15.5 9.5]);
ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

for idx = 1:numel(penalties)
  p = penalties(idx);
  E = Emin + rho.^p .* (E0 - Emin);
  if p == 3
    lineWidth = 2.6;
  else
    lineWidth = 1.25;
  end
  plot(ax, rho, E, 'LineWidth', lineWidth, 'Color', colors(idx, :), ...
    'DisplayName', sprintf('p = %d', p));
end

rhoStar = 0.5;
EStar = Emin + rhoStar^3 * (E0 - Emin);
plot(ax, rhoStar, EStar, 'o', 'MarkerSize', 6.5, ...
  'MarkerFaceColor', colors(3, :), 'MarkerEdgeColor', 'w', ...
  'LineWidth', 0.9, 'HandleVisibility', 'off');
xline(ax, rhoStar, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9, ...
  'HandleVisibility', 'off');
yline(ax, EStar, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9, ...
  'HandleVisibility', 'off');

text(ax, rhoStar + 0.035, EStar + 0.055, ...
  sprintf('V^* = %.1f,  E/E_0 = %.3f', rhoStar, EStar), ...
  'FontSize', 9, 'Color', [0.16 0.16 0.16], ...
  'BackgroundColor', 'w', 'Margin', 2.5);

xlabel(ax, 'Physical density \rho');
ylabel(ax, 'Normalized stiffness E(\rho) / E_0');
title(ax, 'SIMP material interpolation');
legend(ax, 'Location', 'northwest', 'Box', 'off');
xlim(ax, [0 1]);
ylim(ax, [0 1]);
xticks(ax, 0:0.1:1);
yticks(ax, 0:0.1:1);
set(ax, 'FontName', 'Arial', 'FontSize', 9.5, 'LineWidth', 0.8, ...
  'GridAlpha', 0.18, 'MinorGridAlpha', 0.10);

exportgraphics(fig, fullfile(figuresDir, 'simp_interpolation.png'), ...
  'Resolution', 300);
exportgraphics(fig, fullfile(figuresDir, 'simp_interpolation.pdf'), ...
  'ContentType', 'vector');
close(fig);
end
