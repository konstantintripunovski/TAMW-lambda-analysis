%% lambda_3d_figure4_lambda_band_projection.m
% Standalone version of Figure 4 from lambda_3d_visualizations.m.
% Lambda-band section, solid projection column, and projected probability mass.

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
Lambda = lambda_surface(X, Y);

lambda_min = 1 + 3*(ymin/xmax);
lambda_max = 1 + 3*(ymax/xmin);

lambda_demo = 5.5;
c1 = 0.95*lambda_demo;
c2 = 1.05*lambda_demo;

fprintf('lambda_min  = %.6f\n', lambda_min);
fprintf('lambda_max  = %.6f\n', lambda_max);
fprintf('lambda_demo = %.6f\n', lambda_demo);
fprintf('c1          = %.6f\n', c1);
fprintf('c2          = %.6f\n', c2);

%% Global plotting style
set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultTextInterpreter', 'latex');

%% Figure 4: Lambda band, cutting planes, and projected probability mass
fig4 = figure('Name', 'Lambda Band Projection with Probability Density', ...
    'Units', 'pixels', 'Position', [210 90 1120 820]);
ax4 = axes(fig4);
hold(ax4, 'on');

z0 = 0;
plane_alpha = 0.35;
projection_z = z0;
% Display-only lift prevents z-fighting with the probability heatmap.
floor_overlay_z = projection_z + 1e-3;
% Display-only lift separates the highlighted slice from the full surface.
band_surface_offset = 1e-2;
exaggeration_factor = 3;
lambda_cmap = turbo_soft(256);
probability_cmap = hot_modern(256);
lambda_rgb = scalar_to_rgb(Lambda, [lambda_min lambda_max], lambda_cmap);
probability_rgb = scalar_to_rgb(P, [min(P(:)) max(P(:))], probability_cmap);
[X_band, Y_band, Z_band] = lambda_band_patch(c1, c2, xmin, xmax, ymin, ymax, 260, 44);
Lambda_display = z_display(Lambda, lambda_demo, exaggeration_factor);
Z_band_display = z_display(Z_band, lambda_demo, exaggeration_factor);
c1_display = z_display(c1, lambda_demo, exaggeration_factor);
c2_display = z_display(c2, lambda_demo, exaggeration_factor);

% Keep the probability density visible on the floor exactly as in Figure 2.
surf(ax4, X, Y, projection_z*ones(size(P)), probability_rgb, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'texturemap', ...
    'FaceAlpha', 0.96);
surf(ax4, X_band, Y_band, floor_overlay_z*ones(size(Z_band)), ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.15 0.55 1.00], ...
    'FaceAlpha', 0.16);

draw_domain_rectangle(ax4, xmin, xmax, ymin, ymax, z0, [0.08 0.08 0.08], 1.8);

surf(ax4, X, Y, Lambda_display, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.25 0.88 0.98], ...
    'FaceAlpha', 0.35);
draw_lambda_surface_grid(ax4, X, Y, Lambda_display, ...
    [0.20 0.35 0.90], 0.35, 0.35, 13);

plane_x = [xmin xmax; xmin xmax];
plane_y = [ymin ymin; ymax ymax];
plane_c1 = c1_display*ones(size(plane_x));
plane_c2 = c2_display*ones(size(plane_x));
surf(ax4, plane_x, plane_y, plane_c1, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.20 0.47 0.75], ...
    'FaceAlpha', plane_alpha, ...
    'FaceLighting', 'none');
surf(ax4, plane_x, plane_y, plane_c2, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.72 0.28 0.25], ...
    'FaceAlpha', plane_alpha, ...
    'FaceLighting', 'none');

% A solid translucent extrusion makes the band-to-floor projection explicit.
draw_lambda_band_column(ax4, c1, c2, xmin, xmax, ymin, ymax, projection_z, ...
    [0.05 0.40 1.00], 0.15, lambda_demo, exaggeration_factor);

% Continuous blue section matching the reference figure.
surf(ax4, X_band, Y_band, Z_band_display + band_surface_offset, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.15 0.55 1.00], ...
    'FaceAlpha', 0.70);

draw_lambda_level_curve(ax4, c1, xmin, xmax, ymin, ymax, c1_display + band_surface_offset, ...
    'Color', [0.10 0.35 0.75], 'LineWidth', 0.5, 'LineStyle', '-');
draw_lambda_level_curve(ax4, c2, xmin, xmax, ymin, ymax, c2_display + band_surface_offset, ...
    'Color', [0.08 0.30 0.65], 'LineWidth', 0.5, 'LineStyle', '-');
draw_lambda_level_curve(ax4, c1, xmin, xmax, ymin, ymax, floor_overlay_z, ...
    'Color', [0.10 0.35 0.75], 'LineWidth', 0.5, 'LineStyle', '-');
draw_lambda_level_curve(ax4, c2, xmin, xmax, ymin, ymax, floor_overlay_z, ...
    'Color', [0.08 0.30 0.65], 'LineWidth', 0.5, 'LineStyle', '-');

text(ax4, xmin + 0.04*(xmax - xmin), ymax - 0.25*(ymax - ymin), lambda_demo, ...
    sprintf('$\\lambda_{demo}=%.2f$', lambda_demo), ...
    'FontSize', 13, 'Color', [0.18 0.18 0.18]);

colormap(ax4, probability_cmap);
cb4 = colorbar(ax4);
cb4.Label.String = 'Probability density on projected region';
cb4.Label.Interpreter = 'latex';
cb4.FontName = 'Times New Roman';
cb4.TickLabelInterpreter = 'latex';
caxis(ax4, [min(P(:)), max(P(:))]);

xlabel(ax4, '$x$');
ylabel(ax4, '$y$');
zlabel(ax4, '$\lambda$');
title(ax4, 'Lambda-band section and projected probability mass on the $xy$-plane');

view(ax4, 42, 30);
zlim(ax4, [0, 1.18*c2_display]);
format_3d_axes(ax4);
% Keep uniform face colors stable while preserving the heatmap lighting.
set(findobj(ax4, 'Type', 'surface', '-not', 'FaceColor', 'texturemap'), ...
    'FaceLighting', 'none');
annotation(fig4, 'textbox', [0.63 0.83 0.24 0.04], ...
    'String', 'Vertical separation exaggerated for visualization', ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 10, ...
    'Color', [0.28 0.28 0.28], ...
    'EdgeColor', 'none', ...
    'FitBoxToText', 'on');

%% Local helper functions
function P = gaussian_pdf_xy(x, y, mx, sx, my, sy)
P = 1/(2*pi*sx*sy) .* exp(-0.5*((x - mx)./sx).^2 - 0.5*((y - my)./sy).^2);
end

function L = lambda_surface(x, y)
L = 1 + 3*y./x;
end

function Zd = z_display(z_actual, lambda_demo, exaggeration_factor)
Zd = lambda_demo + exaggeration_factor*(z_actual - lambda_demo);
end

function draw_domain_rectangle(ax, xmin, xmax, ymin, ymax, z, color_value, line_width)
plot3(ax, [xmin xmax xmax xmin xmin], ...
    [ymin ymin ymax ymax ymin], ...
    z*ones(1, 5), ...
    'Color', color_value, 'LineWidth', line_width);
end

function draw_lambda_surface_grid(ax, x, y, z, grid_color, grid_alpha, line_width, step_size)
row_idx = unique([1:step_size:size(x, 1), size(x, 1)]);
col_idx = unique([1:step_size:size(x, 2), size(x, 2)]);
mesh(ax, x(row_idx, col_idx), y(row_idx, col_idx), z(row_idx, col_idx), ...
    'FaceColor', 'none', ...
    'EdgeColor', grid_color, ...
    'EdgeAlpha', grid_alpha, ...
    'LineWidth', line_width);
end

function draw_lambda_level_curve(ax, lambda_value, xmin, xmax, ymin, ymax, z, varargin)
slope = (lambda_value - 1)/3;
x_start = max(xmin, ymin/slope);
x_end = min(xmax, ymax/slope);
if x_start < x_end
    x_line = linspace(x_start, x_end, 1200);
    y_line = (x_line/3)*(lambda_value - 1);
    plot3(ax, x_line, y_line, z*ones(size(x_line)), varargin{:});
end
end

function draw_lambda_band_column(ax, c1, c2, xmin, xmax, ymin, ymax, ...
    z_projection, color_value, face_alpha, lambda_demo, exaggeration_factor)
slope1 = (c1 - 1)/3;
slope2 = (c2 - 1)/3;
x_left = max(xmin, ymin/slope2);
x_right = min(xmax, ymax/slope1);

if x_left >= x_right
    return;
end

num_curve = 240;
c1_display = z_display(c1, lambda_demo, exaggeration_factor);
c2_display = z_display(c2, lambda_demo, exaggeration_factor);

% Vertical side walls along the two lambda level-set boundaries.
draw_level_wall(ax, c1, c1_display, xmin, xmax, ymin, ymax, ...
    z_projection, color_value, face_alpha, num_curve);
draw_level_wall(ax, c2, c2_display, xmin, xmax, ymin, ymax, ...
    z_projection, color_value, face_alpha, num_curve);

% If the projected band touches the rectangle, close those domain-edge sides.
draw_domain_edge_wall(ax, 'ymin', c1, c2, xmin, xmax, ymin, ymax, ...
    z_projection, color_value, face_alpha, num_curve, lambda_demo, exaggeration_factor);
draw_domain_edge_wall(ax, 'ymax', c1, c2, xmin, xmax, ymin, ymax, ...
    z_projection, color_value, face_alpha, num_curve, lambda_demo, exaggeration_factor);

% End caps close the extrusion where it meets the rectangular domain.
draw_column_cap(ax, x_left, c1, c2, ymin, ymax, z_projection, ...
    color_value, face_alpha, lambda_demo, exaggeration_factor);
draw_column_cap(ax, x_right, c1, c2, ymin, ymax, z_projection, ...
    color_value, face_alpha, lambda_demo, exaggeration_factor);

end

function [Xb, Yb, Zb] = lambda_band_patch(c1, c2, xmin, xmax, ymin, ymax, num_x, num_eta)
slope1 = (c1 - 1)/3;
slope2 = (c2 - 1)/3;
x_left = max(xmin, ymin/slope2);
x_right = min(xmax, ymax/slope1);

x_values = linspace(x_left, x_right, num_x);
eta = linspace(0, 1, num_eta);
[Xb, Eta] = meshgrid(x_values, eta);

y_low = max(ymin, slope1*x_values);
y_high = min(ymax, slope2*x_values);
Yb = repmat(y_low, num_eta, 1) + Eta.*repmat(y_high - y_low, num_eta, 1);
Zb = lambda_surface(Xb, Yb);
end

function draw_level_wall(ax, lambda_value, lambda_display, xmin, xmax, ...
    ymin, ymax, z_projection, color_value, face_alpha, num_curve)
slope = (lambda_value - 1)/3;
x_start = max(xmin, ymin/slope);
x_end = min(xmax, ymax/slope);

if x_start >= x_end
    return;
end

x_curve = linspace(x_start, x_end, num_curve);
y_curve = slope*x_curve;
surf(ax, [x_curve; x_curve], [y_curve; y_curve], ...
    [z_projection*ones(size(x_curve)); lambda_display*ones(size(x_curve))], ...
    'EdgeColor', 'none', 'FaceColor', color_value, 'FaceAlpha', face_alpha);
end

function draw_domain_edge_wall(ax, edge_name, c1, c2, xmin, xmax, ymin, ...
    ymax, z_projection, color_value, face_alpha, num_curve, lambda_demo, ...
    exaggeration_factor)
slope1 = (c1 - 1)/3;
slope2 = (c2 - 1)/3;

if strcmp(edge_name, 'ymin')
    y_edge = ymin;
    x_start = max(xmin, ymin/slope2);
    x_end = min(xmax, ymin/slope1);
else
    y_edge = ymax;
    x_start = max(xmin, ymax/slope2);
    x_end = min(xmax, ymax/slope1);
end

if x_start >= x_end
    return;
end

x_curve = linspace(x_start, x_end, num_curve);
y_curve = y_edge*ones(size(x_curve));
z_top = z_display(lambda_surface(x_curve, y_curve), lambda_demo, exaggeration_factor);
surf(ax, [x_curve; x_curve], [y_curve; y_curve], ...
    [z_projection*ones(size(x_curve)); z_top], ...
    'EdgeColor', 'none', 'FaceColor', color_value, 'FaceAlpha', face_alpha);
end

function draw_column_cap(ax, x_value, c1, c2, ymin, ymax, z_projection, ...
    color_value, face_alpha, lambda_demo, exaggeration_factor)
slope1 = (c1 - 1)/3;
slope2 = (c2 - 1)/3;
y_low = max(ymin, slope1*x_value);
y_high = min(ymax, slope2*x_value);

if y_low >= y_high
    return;
end

y_values = linspace(y_low, y_high, 90);
t_values = linspace(0, 1, 26);
[Ycap, Tcap] = meshgrid(y_values, t_values);
Xcap = x_value*ones(size(Ycap));
Ztop = z_display(lambda_surface(Xcap, Ycap), lambda_demo, exaggeration_factor);
Zcap = z_projection + Tcap.*(Ztop - z_projection);

surf(ax, Xcap, Ycap, Zcap, ...
    'EdgeColor', 'none', 'FaceColor', color_value, 'FaceAlpha', face_alpha);
end

function rgb = scalar_to_rgb(values, limits, cmap)
t = (values - limits(1))./(limits(2) - limits(1));
t = min(max(t, 0), 1);
idx = 1 + round(t*(size(cmap, 1) - 1));
rgb = reshape(cmap(idx, :), [size(values, 1), size(values, 2), 3]);
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
