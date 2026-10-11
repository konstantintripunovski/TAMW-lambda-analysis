%% lambda_relative_error_80.m
% Relative-error interval around the peak lambda giving 80% probability.
%
% Implementation plan:
% 1) Define the bounded physiological region R and the independent Gaussian
%    density P(x,y).
% 2) Compute MR = integral_R P(x,y) dx dy numerically with integral2.
% 3) Compute P_delta(lambda) inside the original band-probability framework:
%    for each x, intersect the admissible y-band with [ymin,ymax], integrate
%    the Gaussian y-density over that band, then integrate over x.
% 4) Evaluate P_0.05(lambda) over a fine lambda grid on the full admissible
%    lambda support.
% 5) Find lambda_peak from the grid maximum and refine it with fminbnd.
% 6) With lambda fixed at lambda_peak, compute Q(delta) = P_delta(lambda_peak)
%    over a fine delta grid and solve Q(delta) = 0.80 with fzero.
% 7) Compute and display [lambda_low, lambda_high], then create the three
%    requested manuscript-friendly figures.

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

delta0 = 0.05;
target_probability = 0.80;

num_lambda_grid = 3000;
num_delta_grid = 2500;
integration_abs_tol = 1e-11;
integration_rel_tol = 1e-9;

%% Lambda support and normalization over R
lambda_min = 1 + 3*(ymin/xmax);
lambda_max = 1 + 3*(ymax/xmin);

MR = compute_MR(mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

assert(MR > 0 && MR < 1, 'MR should be between 0 and 1.');

%% P_0.05(lambda) over the full admissible lambda range
lambda_grid = linspace(lambda_min, lambda_max, num_lambda_grid);
P_delta0_grid = compute_Pcurve_for_fixed_delta(lambda_grid, delta0, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

check_numeric_array(P_delta0_grid, 'P_0.05(lambda) grid');
check_probability_bounds(P_delta0_grid, 'P_0.05(lambda) grid');

[P_peak_grid, idx_peak] = max(P_delta0_grid);
lambda_peak_grid = lambda_grid(idx_peak);

% Refine the peak near the grid maximum while keeping the grid-based curve
% as the primary object shown in Figure 1.
idx_left = max(idx_peak - 2, 1);
idx_right = min(idx_peak + 2, numel(lambda_grid));
peak_objective = @(lambda_value) -compute_Pdelta(lambda_value, delta0, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

if idx_left < idx_right
    lambda_peak = fminbnd(peak_objective, lambda_grid(idx_left), lambda_grid(idx_right));
    P_peak = -peak_objective(lambda_peak);
else
    lambda_peak = lambda_peak_grid;
    P_peak = P_peak_grid;
end

%% Q(delta) = P_delta(lambda_peak)
delta_full_support = max(1 - lambda_min/lambda_peak, lambda_max/lambda_peak - 1);
delta_max = max(1.0, 1.05*delta_full_support);
delta_grid = linspace(0, delta_max, num_delta_grid);

Q_grid = compute_Qcurve_for_fixed_lambda(delta_grid, lambda_peak, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

check_numeric_array(Q_grid, 'Q(delta) grid');
check_probability_bounds(Q_grid, 'Q(delta) grid');

delta_80 = find_delta80(delta_grid, Q_grid, target_probability, lambda_peak, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

P_check = compute_Pdelta(lambda_peak, delta_80, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, ...
    integration_abs_tol, integration_rel_tol);

lambda_low = (1 - delta_80)*lambda_peak;
lambda_high = (1 + delta_80)*lambda_peak;

if abs(P_check - target_probability) > 5e-5
    warning('Final probability check differs from %.4f by %.4g.', ...
        target_probability, abs(P_check - target_probability));
end

%% Print results
fprintf('\n============================================================\n');
fprintf('Relative-error interval around lambda_peak for 80%% probability\n');
fprintf('============================================================\n');
fprintf('MR                                  = %.12f\n', MR);
fprintf('lambda_min                          = %.12f\n', lambda_min);
fprintf('lambda_max                          = %.12f\n', lambda_max);
fprintf('lambda_peak                         = %.12f\n', lambda_peak);
fprintf('P_0.05(lambda_peak)                 = %.12f\n', P_peak);
fprintf('delta_80                            = %.12f\n', delta_80);
fprintf('lambda_low                          = %.12f\n', lambda_low);
fprintf('lambda_high                         = %.12f\n', lambda_high);
fprintf('Final interval [lambda_low,lambda_high] = [%.12f, %.12f]\n', ...
    lambda_low, lambda_high);
fprintf('Numerical check P_delta_80(lambda_peak) = %.12f\n', P_check);
fprintf('Target probability                  = %.12f\n', target_probability);
fprintf('Absolute check error                = %.3e\n', abs(P_check - target_probability));
fprintf('============================================================\n\n');

%% Figure styling
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 13);
set(groot, 'defaultLineLineWidth', 2.0);

%% Figure 1: P_0.05(lambda), peak, and final interval
figure(1);
hold on;
grid on;
box on;

interval_color = [0.80 0.90 1.00];
curve_color = [0.00 0.25 0.55];
peak_color = [0.75 0.10 0.10];
edge_color = [0.10 0.45 0.75];

y_shade_top = max(P_delta0_grid)*1.08;
patch([lambda_low lambda_high lambda_high lambda_low], ...
      [0 0 y_shade_top y_shade_top], interval_color, ...
      'EdgeColor', 'none', 'FaceAlpha', 0.45);

plot(lambda_grid, P_delta0_grid, '-', 'Color', curve_color, 'LineWidth', 2.3);
plot(lambda_peak, P_peak, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', peak_color, 'MarkerEdgeColor', peak_color);
xline(lambda_peak, '-', 'Color', peak_color, 'LineWidth', 1.8);
xline(lambda_low, '--', 'Color', edge_color, 'LineWidth', 1.8);
xline(lambda_high, '--', 'Color', edge_color, 'LineWidth', 1.8);

xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('$P_{0.05}(\lambda)$', 'Interpreter', 'latex');
title('$P_{0.05}(\lambda)$ and the 80\% interval around $\lambda_{\rm peak}$', ...
    'Interpreter', 'latex');
ylim([0 y_shade_top]);
xlim([lambda_min 15]);

label_y = y_shade_top*0.055;
text(lambda_low, label_y, sprintf('$\\lambda_{low}=%.4f$', lambda_low), ...
    'Interpreter', 'latex', 'Rotation', 90, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);
text(lambda_peak, label_y, sprintf('$\\lambda_{peak}=%.4f$', lambda_peak), ...
    'Interpreter', 'latex', 'Rotation', 90, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);
text(lambda_high, label_y, sprintf('$\\lambda_{high}=%.4f$', lambda_high), ...
    'Interpreter', 'latex', 'Rotation', 90, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12);

annotation_text = sprintf('$80\\%%$ interval around $\\lambda_{\\mathrm{peak}}$: [%.4f, %.4f]', ...
    lambda_low, lambda_high);
text(lambda_min + 0.03*(lambda_max - lambda_min), y_shade_top*0.92, annotation_text, ...
    'Interpreter', 'latex', 'FontSize', 12, ...
    'BackgroundColor', 'w', 'Margin', 6, 'EdgeColor', [0.75 0.75 0.75]);

hold off;

%% Figure 2: Q(delta) and delta_80
figure(2);
hold on;
grid on;
box on;

plot(delta_grid, Q_grid, '-', 'Color', curve_color, 'LineWidth', 2.3);
yline(target_probability, '--', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.7);
xline(delta_80, '--', 'Color', peak_color, 'LineWidth', 1.7);
plot(delta_80, target_probability, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', peak_color, 'MarkerEdgeColor', peak_color);

xlabel('$\delta$', 'Interpreter', 'latex');
ylabel('$Q(\delta)=P_{\delta}(\lambda_{\rm peak})$', 'Interpreter', 'latex');
title('Finding the relative error $\delta_{80}$ from $Q(\delta)$', ...
    'Interpreter', 'latex');
ylim([0 1.03]);
xlim([0 min(delta_max, max(delta_80*1.55, delta_80 + 0.05))]);

text(delta_80, target_probability + 0.045, ...
    sprintf('$(\\delta_{80},0.80)=(%.4f,0.80)$', delta_80), ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, ...
    'BackgroundColor', 'w', 'Margin', 5);

hold off;

%% Figure 3: Minimal final interval figure
figure(3);
hold on;
grid on;
box on;

plot([lambda_low lambda_high], [1 1], '-', 'Color', curve_color, 'LineWidth', 7);
plot(lambda_peak, 1, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', peak_color, 'MarkerEdgeColor', peak_color);
plot([lambda_low lambda_low], [0.92 1.08], '-', 'Color', edge_color, 'LineWidth', 2);
plot([lambda_high lambda_high], [0.92 1.08], '-', 'Color', edge_color, 'LineWidth', 2);

text(lambda_low, 1.12, sprintf('$\\lambda_{low}=%.4f$', lambda_low), ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(lambda_peak, 0.84, sprintf('$\\lambda_{peak}=%.4f$', lambda_peak), ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(lambda_high, 1.12, sprintf('$\\lambda_{high}=%.4f$', lambda_high), ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);

xlabel('$\lambda$', 'Interpreter', 'latex');
yticks([]);
ylim([0.70 1.30]);
padding = 0.12*(lambda_high - lambda_low);
xlim([lambda_low - padding, lambda_high + padding]);
title('Resulting lambda interval containing 80\% of the population around $\lambda_{\rm peak}$', ...
    'Interpreter', 'latex');

hold off;

%% Local helper functions
function p = gaussian_pdf_xy(x, y, mx, sx, my, sy)
    p = 1./(2*pi*sx*sy) .* exp(-0.5*((x - mx)./sx).^2 ...
        - 0.5*((y - my)./sy).^2);
end

function gx = gaussian_pdf_1d(x, mu, sigma)
    gx = 1./(sqrt(2*pi)*sigma) .* exp(-0.5*((x - mu)./sigma).^2);
end

function F = gaussian_cdf_1d(z, mu, sigma)
    F = 0.5 .* (1 + erf((z - mu)./(sigma*sqrt(2))));
end

function MR = compute_MR(mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol)
    pdf_handle = @(x,y) gaussian_pdf_xy(x, y, mx, sx, my, sy);
    MR = integral2(pdf_handle, xmin, xmax, ymin, ymax, ...
        'AbsTol', abs_tol, 'RelTol', rel_tol, 'Method', 'iterated');
end

function Pdelta = compute_Pdelta(lambda_value, delta, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    integrand = @(x) x_integrand_for_band(x, lambda_value, delta, ...
        mx, sx, my, sy, ymin, ymax);

    raw_probability = integral(integrand, xmin, xmax, ...
        'ArrayValued', true, 'AbsTol', abs_tol, 'RelTol', rel_tol);

    Pdelta = raw_probability ./ MR;
    Pdelta = min(max(Pdelta, 0), 1);
end

function values = x_integrand_for_band(x, lambda_value, delta, ...
    mx, sx, my, sy, ymin, ymax)

    band_lower = (x./3) .* ((1 - delta).*lambda_value - 1);
    band_upper = (x./3) .* ((1 + delta).*lambda_value - 1);

    y_lower = max(ymin, band_lower);
    y_upper = min(ymax, band_upper);

    valid_band = y_lower < y_upper;
    y_probability = zeros(size(x));

    if any(valid_band(:))
        y_probability(valid_band) = ...
            gaussian_cdf_1d(y_upper(valid_band), my, sy) ...
            - gaussian_cdf_1d(y_lower(valid_band), my, sy);
    end

    values = gaussian_pdf_1d(x, mx, sx) .* y_probability;
end

function Pcurve = compute_Pcurve_for_fixed_delta(lambda_values, delta, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    Pcurve = arrayfun(@(lambda_value) compute_Pdelta(lambda_value, delta, MR, ...
        mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol), lambda_values);
end

function Qcurve = compute_Qcurve_for_fixed_lambda(delta_values, lambda_value, MR, ...
    mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    Qcurve = arrayfun(@(delta) compute_Pdelta(lambda_value, delta, MR, ...
        mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol), delta_values);
end

function delta_80 = find_delta80(delta_grid, Q_grid, target_probability, ...
    lambda_peak, MR, mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    if Q_grid(1) > target_probability
        error('Q(0) already exceeds the target probability.');
    end

    if Q_grid(end) < target_probability
        error('The delta grid does not reach the target probability.');
    end

    first_idx = find(Q_grid >= target_probability, 1, 'first');

    if first_idx == 1
        bracket = [delta_grid(1), delta_grid(2)];
    else
        bracket = [delta_grid(first_idx - 1), delta_grid(first_idx)];
    end

    root_function = @(delta) compute_Pdelta(lambda_peak, delta, MR, ...
        mx, sx, my, sy, xmin, xmax, ymin, ymax, abs_tol, rel_tol) ...
        - target_probability;

    delta_80 = fzero(root_function, bracket);
end

function check_numeric_array(values, name)
    if any(isnan(values(:))) || any(isinf(values(:)))
        error('%s contains NaN or Inf values.', name);
    end
end

function check_probability_bounds(values, name)
    tol = 1e-8;
    if any(values(:) < -tol) || any(values(:) > 1 + tol)
        error('%s contains values outside [0,1] beyond tolerance.', name);
    end
end
