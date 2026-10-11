%% lambda_3d_figure1_gaussian_marginals.m
% Standalone version of Figure 1 from lambda_3d_visualizations.m.
% Joint Gaussian surface with the x- and y-marginal Gaussians.

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
Px = gaussian_pdf_x(x, mx, sx);
Py = gaussian_pdf_y(y, my, sy);

%% Global plotting style
set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultTextInterpreter', 'latex');

%% Figure 1: Joint Gaussian surface with marginal Gaussians
fig1 = figure('Name', 'Joint Gaussian with Marginals', ...
    'Units', 'pixels', 'Position', [80 80 1050 780]);
ax1 = axes(fig1);
hold(ax1, 'on');

surf(ax1, X, Y, P, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceAlpha', 0.92);

z_floor_1 = 0;
marginal_scale = 0.85*max(P(:))/max([Px(:); Py(:)]);

plot3(ax1, x, ymin*ones(size(x)), z_floor_1 + marginal_scale*Px, ...
    'Color', [0.05 0.20 0.55], 'LineWidth', 2.2);
plot3(ax1, xmin*ones(size(y)), y, z_floor_1 + marginal_scale*Py, ...
    'Color', [0.58 0.12 0.08], 'LineWidth', 2.2);
plot3(ax1, x, ymin*ones(size(x)), z_floor_1*ones(size(x)), ...
    'Color', [0.05 0.20 0.55], 'LineWidth', 1.4, 'LineStyle', ':');
plot3(ax1, xmin*ones(size(y)), y, z_floor_1*ones(size(y)), ...
    'Color', [0.58 0.12 0.08], 'LineWidth', 1.4, 'LineStyle', ':');

draw_domain_rectangle(ax1, xmin, xmax, ymin, ymax, z_floor_1, [0.18 0.18 0.18], 1.6);

colormap(ax1, parula(256));
cb1 = colorbar(ax1);
cb1.Label.String = '$P(x,y)$';
cb1.Label.Interpreter = 'latex';
cb1.FontName = 'Times New Roman';
cb1.TickLabelInterpreter = 'latex';

xlabel(ax1, '$x$');
ylabel(ax1, '$y$');
zlabel(ax1, 'Probability density');
title(ax1, 'Joint Gaussian density with its $x$- and $y$-marginals');
legend(ax1, {'$P(x,y)$', '$P_x(x)$ along $y=y_{\min}$', ...
    '$P_y(y)$ along $x=x_{\min}$'}, ...
    'Location', 'northeast', 'Box', 'off');

view(ax1, 42, 28);
zlim(ax1, [z_floor_1, 1.08*max(P(:))]);
format_3d_axes(ax1);

%% Local helper functions
function P = gaussian_pdf_xy(x, y, mx, sx, my, sy)
P = 1/(2*pi*sx*sy) .* exp(-0.5*((x - mx)./sx).^2 - 0.5*((y - my)./sy).^2);
end

function P = gaussian_pdf_x(x, mx, sx)
P = 1/(sqrt(2*pi)*sx) .* exp(-0.5*((x - mx)./sx).^2);
end

function P = gaussian_pdf_y(y, my, sy)
P = 1/(sqrt(2*pi)*sy) .* exp(-0.5*((y - my)./sy).^2);
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
axes(ax); %#ok<LAXES>
lighting gouraud;
camlight headlight;
material dull;
end

function set_concise_xy_ticks(ax)
xticks(ax, [50 150]);
yticks(ax, [50 150 250]);
end
