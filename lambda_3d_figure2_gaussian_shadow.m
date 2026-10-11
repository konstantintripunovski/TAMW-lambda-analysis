%% lambda_3d_figure2_gaussian_shadow.m
% Standalone version of Figure 2 from lambda_3d_visualizations.m.
% Probability-density shadow of the joint Gaussian on the xy-plane.

clear; clc; close all;

%% Parameters
mx = 97;
sx = 30;
my = 139;
sy = 31;

xmin = 19;
xmax = 195;
ymin = 57;
ymax = 245;

num_x = 260;
num_y = 260;

x = linspace(xmin, xmax, num_x);
y = linspace(ymin, ymax, num_y);
[X, Y] = meshgrid(x, y);

P = gaussian_pdf_xy(X, Y, mx, sx, my, sy);

%% Global plotting style
set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultTextInterpreter', 'latex');

%% Figure 2: Probability-density shadow on the xy-plane
fig2 = figure('Name', 'Gaussian Shadow on xy-plane', ...
    'Units', 'pixels', 'Position', [130 100 980 760]);
ax2 = axes(fig2);
hold(ax2, 'on');

surf(ax2, X, Y, zeros(size(P)), P, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.98);
draw_domain_rectangle(ax2, xmin, xmax, ymin, ymax, 0, [0.08 0.08 0.08], 2.0);

colormap(ax2, hot_modern(256));
cb2 = colorbar(ax2);
cb2.Label.String = 'Probability density';
cb2.Label.Interpreter = 'latex';
cb2.FontName = 'Times New Roman';
cb2.TickLabelInterpreter = 'latex';

xlabel(ax2, '$x$');
ylabel(ax2, '$y$');
zlabel(ax2, '$z=0$');
title(ax2, 'Projection of the joint Gaussian density onto the $xy$-plane');

view(ax2, 37, 36);
zlim(ax2, [0, 1]);
zticks(ax2, 0);
format_3d_axes(ax2);

%% Local helper functions
function P = gaussian_pdf_xy(x, y, mx, sx, my, sy)
P = 1/(2*pi*sx*sy) .* exp(-0.5*((x - mx)./sx).^2 - 0.5*((y - my)./sy).^2);
end

function draw_domain_rectangle(ax, xmin, xmax, ymin, ymax, z, color_value, line_width)
plot3(ax, [xmin xmax xmax xmin xmin], ...
    [ymin ymin ymax ymax ymin], ...
    z*ones(1, 5), ...
    'Color', color_value, 'LineWidth', line_width);
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
axes(ax); 
lighting gouraud;
camlight headlight;
material dull;
end

function set_concise_xy_ticks(ax)
xticks(ax, [50 150]);
yticks(ax, [50 150 250]);
end

function cmap = hot_modern(n)
if nargin < 1
    n = 256;
end
base = [ ...
    0.98 0.98 0.96
    0.82 0.92 0.86
    0.50 0.76 0.70
    0.28 0.56 0.68
    0.34 0.22 0.48];
t = linspace(0, 1, size(base, 1));
tq = linspace(0, 1, n);
cmap = interp1(t, base, tq, 'pchip');
cmap = min(max(cmap, 0), 1);
end
