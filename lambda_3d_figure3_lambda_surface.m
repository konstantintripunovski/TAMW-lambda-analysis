%% lambda_3d_figure3_lambda_surface.m
% Standalone version of Figure 3 from lambda_3d_visualizations.m.
% Curved lambda surface over the admissible xy-domain.

clear; clc; close all;

%% Parameters
xmin = 19;
xmax = 195;
ymin = 57;
ymax = 245;

num_x = 260;
num_y = 260;

x = linspace(xmin, xmax, num_x);
y = linspace(ymin, ymax, num_y);
[X, Y] = meshgrid(x, y);

Lambda = lambda_surface(X, Y);

lambda_min = 1 + 3*(ymin/xmax);
lambda_max = 1 + 3*(ymax/xmin);

%% Global plotting style
set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultTextInterpreter', 'latex');

%% Figure 3: Curved lambda surface
fig3 = figure('Name', 'Lambda Surface', ...
    'Units', 'pixels', 'Position', [180 120 1030 780]);
ax3 = axes(fig3);
hold(ax3, 'on');

surf(ax3, X, Y, Lambda, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.94);

plot3(ax3, xmax, ymin, lambda_min, 'o', ...
    'MarkerFaceColor', [0.02 0.35 0.45], ...
    'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4);
plot3(ax3, xmin, ymax, lambda_max, 'o', ...
    'MarkerFaceColor', [0.65 0.08 0.08], ...
    'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4);
text(ax3, xmax, ymin, lambda_min, ...
    sprintf('  $\\lambda_{\\min}=%.2f$', lambda_min), ...
    'FontSize', 13, 'Color', [0.02 0.22 0.30]);
text(ax3, xmin, ymax, lambda_max, ...
    sprintf('  $\\lambda_{\\max}=%.2f$', lambda_max), ...
    'FontSize', 13, 'Color', [0.42 0.04 0.04]);

colormap(ax3, turbo_soft(256));
cb3 = colorbar(ax3);
cb3.Label.String = '$\lambda(x,y)$';
cb3.Label.Interpreter = 'latex';
cb3.FontName = 'Times New Roman';
cb3.TickLabelInterpreter = 'latex';

xlabel(ax3, '$x$');
ylabel(ax3, '$y$');
zlabel(ax3, '$\lambda$');
title(ax3, 'Curved lambda surface over the admissible $xy$-domain');

view(ax3, 42, 30);
zlim(ax3, [lambda_min, lambda_max]);
format_3d_axes(ax3);

%% Local helper functions
function L = lambda_surface(x, y)
L = 1 + 3*y./x;
end

function format_3d_axes(ax)
set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 13, ...
    'LineWidth', 1.25, ...
    'TickDir', 'out', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'ZMinorTick', 'on', ...
    'GridAlpha', 0.18, ...
    'MinorGridAlpha', 0.08);
grid(ax, 'on');
set_concise_xy_ticks(ax);
pbaspect(ax, [1.25 1 0.72]);
axes(ax); %#ok<LAXES>
lighting gouraud;
camlight headlight;
material dull;
end

function set_concise_xy_ticks(ax)
xticks(ax, [50 150]);
yticks(ax, [50 150 250]);
end

function cmap = turbo_soft(n)
if nargin < 1
    n = 256;
end
try
    cmap = turbo(n);
catch
    cmap = parula(n);
end
cmap = 0.88*cmap + 0.12;
cmap = min(max(cmap, 0), 1);
end
