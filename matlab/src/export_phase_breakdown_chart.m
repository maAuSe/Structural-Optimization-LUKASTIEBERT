function export_phase_breakdown_chart()
%EXPORT_PHASE_BREAKDOWN_CHART Create the report figure for timing phases.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
figuresDir = fullfile(repoRoot, 'figures');
if ~isfolder(figuresDir)
  mkdir(figuresDir);
end

caseNames = {
  'Classical OC sens.'
  'MTOP OC sens.'
  'Classical OC dens.'
  'MTOP OC dens.'
  'MTOP MMA sens.'
  'MTOP MMA Heav. \eta=0.5'
};

totalTime = [263.23; 12.28; 692.91; 136.57; 187.28; 370.87];
phaseTime = [
  43.03   214.77  1.02   0.91
  4.33    5.35    1.01   0.88
  102.05  498.31  7.32   76.88
  12.16   16.03   9.17   96.84
  10.95   13.45   2.79   158.01
  29.13   36.72   23.13  276.62
];
residualTime = max(totalTime - sum(phaseTime, 2), 0);
plotTime = [phaseTime residualTime];

phaseNames = {'FE asm.', 'FE solve', 'Filter / chain rule', ...
  'Optimizer', 'Other'};
phaseColors = [
  0.24 0.45 0.70
  0.12 0.62 0.56
  0.93 0.68 0.19
  0.79 0.29 0.23
  0.70 0.70 0.70
];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
  'Position', [2 2 18.2 12.4]);
layout = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', ...
  'Padding', 'compact');
ax = nexttile(layout);
hold(ax, 'on');
box(ax, 'off');
grid(ax, 'on');

bars = barh(ax, plotTime, 0.68, 'stacked', 'EdgeColor', 'none');
for idx = 1:numel(bars)
  bars(idx).FaceColor = phaseColors(idx, :);
end

set(ax, 'YDir', 'reverse', 'YTick', 1:numel(caseNames), ...
  'YTickLabel', caseNames, 'FontName', 'Arial', 'FontSize', 8.8, ...
  'LineWidth', 0.8, 'GridAlpha', 0.16, 'MinorGridAlpha', 0.08, ...
  'TickDir', 'out');
xlabel(ax, 'Cumulative wall-clock time [s]');
title(ax, 'Runtime composition by algorithmic phase', ...
  'FontWeight', 'bold', 'FontSize', 11);
xlim(ax, [0 735]);
xticks(ax, 0:100:700);

legend(ax, bars, phaseNames, 'Location', 'southoutside', ...
  'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8.2);

for row = 1:numel(totalTime)
  text(ax, totalTime(row) + 9, row, sprintf('%.1f s', totalTime(row)), ...
    'FontName', 'Arial', 'FontSize', 8.1, 'Color', [0.12 0.12 0.12], ...
    'VerticalAlignment', 'middle');
end

% Label the dominant phase in each row where the segment is wide enough.
for row = 1:size(plotTime, 1)
  cumulative = 0;
  for col = 1:size(plotTime, 2)
    width = plotTime(row, col);
    if width > 44
      text(ax, cumulative + 0.5 * width, row, sprintf('%.0f%%', ...
        100 * width / totalTime(row)), 'FontName', 'Arial', ...
        'FontSize', 7.8, 'FontWeight', 'bold', 'Color', 'w', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    cumulative = cumulative + width;
  end
end

text(ax, 372, 1.50, 'MTOP removes most FE-solve time', ...
  'FontName', 'Arial', 'FontSize', 8.3, 'FontWeight', 'bold', ...
  'Color', [0.05 0.32 0.30], 'BackgroundColor', 'w', 'Margin', 2.5);
plot(ax, [214.77 5.35], [1 2], '-', 'Color', [0.05 0.32 0.30], ...
  'LineWidth', 1.2, 'HandleVisibility', 'off');
text(ax, 112, 1.54, '40x', 'FontName', 'Arial', 'FontSize', 8.1, ...
  'Color', [0.05 0.32 0.30], 'BackgroundColor', 'w', 'Margin', 1.8);

text(ax, 410, 5.52, 'MMA is optimizer-dominated', ...
  'FontName', 'Arial', 'FontSize', 8.3, 'FontWeight', 'bold', ...
  'Color', [0.50 0.13 0.10], 'BackgroundColor', 'w', 'Margin', 2.5);
plot(ax, [158.01 276.62], [5 6], '-', 'Color', [0.50 0.13 0.10], ...
  'LineWidth', 1.2, 'HandleVisibility', 'off');

ylim(ax, [0.35 numel(caseNames) + 0.65]);

exportgraphics(fig, fullfile(figuresDir, 'phase_breakdown_chart.png'), ...
  'Resolution', 300);
exportgraphics(fig, fullfile(figuresDir, 'phase_breakdown_chart.pdf'), ...
  'ContentType', 'vector');
close(fig);
end
