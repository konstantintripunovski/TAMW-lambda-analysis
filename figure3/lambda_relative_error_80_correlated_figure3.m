%% lambda_relative_error_80.m
% Relative-error interval around the peak lambda giving 80% probability.
%
% Implementation plan:
% 1) Define the bounded physiological region R and the correlated bivariate
%    Gaussian density P_rho(x,y).
% 2) For each rho, compute MR(rho) = integral_R P_rho(x,y) dx dy.
% 3) Compute P_delta(lambda) inside the original band-probability framework:
%    for each x, intersect the admissible y-band with [ymin,ymax], integrate
%    the conditional Gaussian Y|X=x over that band, then integrate over x.
% 4) Evaluate P_0.05(lambda) over a fine lambda grid on the full admissible
%    lambda support.
% 5) Find lambda_peak from the grid maximum and refine it with fminbnd.
% 6) With lambda fixed at lambda_peak, compute Q(delta) = P_delta(lambda_peak)
%    over a fine delta grid and solve Q(delta) = 0.80 with fzero.
% 7) Compute and display [lambda_low, lambda_high], then create only the
%    final interval comparison figure for the requested rho values.

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

rho_values = [0, 0.25, 0.50, 0.75];

delta0 = 0.05;
target_probability = 0.80;

num_lambda_grid = 3000;
num_delta_grid = 2500;
integration_abs_tol = 1e-11;
integration_rel_tol = 1e-9;

%% Lambda support
lambda_min = 1 + 3*(ymin/xmax);
lambda_max = 1 + 3*(ymax/xmin);
lambda_grid = linspace(lambda_min, lambda_max, num_lambda_grid);

%% Complete analysis for each rho
num_rho = numel(rho_values);
results = repmat(struct( ...
    'rho', NaN, ...
    'MR', NaN, ...
    'lambda_peak', NaN, ...
    'P_peak', NaN, ...
    'delta_80', NaN, ...
    'lambda_low', NaN, ...
    'lambda_high', NaN, ...
    'P_check', NaN, ...
    'delta_grid', [], ...
    'Q_grid', [], ...
    'P_delta0_grid', []), num_rho, 1);

for k = 1:num_rho
    rho = rho_values(k);

    MR = compute_MR(mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    assert(MR > 0 && MR < 1, 'MR should be between 0 and 1 for rho = %.2f.', rho);

    P_delta0_grid = compute_Pcurve_for_fixed_delta(lambda_grid, delta0, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    check_numeric_array(P_delta0_grid, sprintf('P_0.05(lambda) grid, rho = %.2f', rho));
    check_probability_bounds(P_delta0_grid, sprintf('P_0.05(lambda) grid, rho = %.2f', rho));

    [P_peak_grid, idx_peak] = max(P_delta0_grid);
    lambda_peak_grid = lambda_grid(idx_peak);

    % Refine the peak near the grid maximum while keeping the grid-based
    % curve available for the same analysis as the full correlated script.
    idx_left = max(idx_peak - 2, 1);
    idx_right = min(idx_peak + 2, numel(lambda_grid));
    peak_objective = @(lambda_value) -compute_Pdelta(lambda_value, delta0, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    if idx_left < idx_right
        lambda_peak = fminbnd(peak_objective, lambda_grid(idx_left), lambda_grid(idx_right));
        P_peak = -peak_objective(lambda_peak);
    else
        lambda_peak = lambda_peak_grid;
        P_peak = P_peak_grid;
    end

    delta_full_support = max(1 - lambda_min/lambda_peak, lambda_max/lambda_peak - 1);
    delta_max = max(1.0, 1.05*delta_full_support);
    delta_grid = linspace(0, delta_max, num_delta_grid);

    Q_grid = compute_Qcurve_for_fixed_lambda(delta_grid, lambda_peak, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    check_numeric_array(Q_grid, sprintf('Q(delta) grid, rho = %.2f', rho));
    check_probability_bounds(Q_grid, sprintf('Q(delta) grid, rho = %.2f', rho));

    delta_80 = find_delta80(delta_grid, Q_grid, target_probability, lambda_peak, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    P_check = compute_Pdelta(lambda_peak, delta_80, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, ...
        integration_abs_tol, integration_rel_tol);

    lambda_low = (1 - delta_80)*lambda_peak;
    lambda_high = (1 + delta_80)*lambda_peak;

    if abs(P_check - target_probability) > 5e-5
        warning('For rho = %.2f, final probability check differs from %.4f by %.4g.', ...
            rho, target_probability, abs(P_check - target_probability));
    end

    results(k).rho = rho;
    results(k).MR = MR;
    results(k).lambda_peak = lambda_peak;
    results(k).P_peak = P_peak;
    results(k).delta_80 = delta_80;
    results(k).lambda_low = lambda_low;
    results(k).lambda_high = lambda_high;
    results(k).P_check = P_check;
    results(k).delta_grid = delta_grid;
    results(k).Q_grid = Q_grid;
    results(k).P_delta0_grid = P_delta0_grid;
end

%% Print results
fprintf('\n============================================================\n');
fprintf('Relative-error interval around lambda_peak for 80%% probability\n');
fprintf('Correlated bivariate normal analysis\n');
fprintf('============================================================\n');
fprintf('lambda_min                          = %.12f\n', lambda_min);
fprintf('lambda_max                          = %.12f\n', lambda_max);
fprintf('Target probability                  = %.12f\n', target_probability);
fprintf('------------------------------------------------------------\n');

for k = 1:num_rho
    fprintf('rho                                 = %.2f\n', results(k).rho);
    fprintf('MR                                  = %.12f\n', results(k).MR);
    fprintf('lambda_peak                         = %.12f\n', results(k).lambda_peak);
    fprintf('P_0.05(lambda_peak)                 = %.12f\n', results(k).P_peak);
    fprintf('delta_80                            = %.12f\n', results(k).delta_80);
    fprintf('lambda_low                          = %.12f\n', results(k).lambda_low);
    fprintf('lambda_high                         = %.12f\n', results(k).lambda_high);
    fprintf('Final interval [lambda_low,lambda_high] = [%.12f, %.12f]\n', ...
        results(k).lambda_low, results(k).lambda_high);
    fprintf('Numerical check P_delta_80(lambda_peak) = %.12f\n', results(k).P_check);
    fprintf('Absolute check error                = %.3e\n', ...
        abs(results(k).P_check - target_probability));
    fprintf('------------------------------------------------------------\n');
end

fprintf('============================================================\n\n');

%% Figure styling
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 13);
set(groot, 'defaultLineLineWidth', 2.0);

curve_colors = [ ...
    0.00 0.25 0.55
    0.00 0.50 0.45
    0.78 0.36 0.07
    0.62 0.12 0.22];
peak_color = [0.75 0.10 0.10];
edge_color = [0.10 0.45 0.75];

%% Figure 3: Final 80% interval comparison
figure(3);
clf;
hold on;

% Put rho = 0 at the top
y_positions = num_rho:-1:1;

% Optional reference line through the rho = 0 peak
xline(results(1).lambda_peak, ':', ...
    'Color', [0.70 0.70 0.70], ...
    'LineWidth', 1.0, ...
    'HandleVisibility', 'off');

% Formatting controls
cap_halfheight = 0.08;     % shorter endpoint caps
cap_linewidth = 1.3;       % thinner endpoint caps
label_offset = 0.12;       % height of numerical labels above bar
rho_offset = 0.31;         % height of rho labels above bar
rho_x = 4.15;              % same horizontal location for every rho label

for k = 1:num_rho
    y0 = y_positions(k);

    % Main 80% interval
    plot([results(k).lambda_low results(k).lambda_high], ...
        [y0 y0], '-', ...
        'Color', curve_colors(k, :), ...
        'LineWidth', 7);

    % Peak
    plot(results(k).lambda_peak, y0, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', peak_color, 'MarkerEdgeColor', peak_color);

    % Endpoint caps
    plot([results(k).lambda_low results(k).lambda_low], ...
        [y0 - cap_halfheight y0 + cap_halfheight], '-', ...
        'Color', edge_color, ...
        'LineWidth', cap_linewidth);

    plot([results(k).lambda_high results(k).lambda_high], ...
        [y0 - cap_halfheight y0 + cap_halfheight], '-', ...
        'Color', edge_color, ...
        'LineWidth', cap_linewidth);

    % Numerical labels: 3 decimal places, no "X"
    text(results(k).lambda_low, y0 + label_offset, ...
        sprintf('%.3f', results(k).lambda_low), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.70 0.70 0.70], ...
        'Margin', 3);

    text(results(k).lambda_peak, y0 + label_offset, ...
        sprintf('%.3f', results(k).lambda_peak), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.70 0.70 0.70], ...
        'Margin', 3);

    text(results(k).lambda_high, y0 + label_offset, ...
        sprintf('%.3f', results(k).lambda_high), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.70 0.70 0.70], ...
        'Margin', 3);

    % rho labels, all aligned at the same x-position
    text(rho_x, y0 + rho_offset, ...
        sprintf('$\\rho = %.2g$', results(k).rho), ...
        'Interpreter', 'latex', ...
        'FontSize', 15, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end

%% Axes
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 15);

% rho is already written directly on the graph
yticks([]);

ylim([0.5 num_rho + 0.55]);
xlim([2.5 7.5]);

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 13;

% Only vertical grid lines
ax.XGrid = 'on';
ax.YGrid = 'off';

% Lighter grid
ax.GridColor = [0.75 0.75 0.75];
ax.GridAlpha = 0.25;

% Removes the mirrored upper/right ticks and frame
box off;

% No title -- caption can explain the figure

hold off;

%% Local helper functions
function p = gaussian_pdf_xy(x, y, mx, sx, my, sy, rho)
    one_minus_rho2 = 1 - rho^2;
    z = ((x - mx)./sx).^2 ...
        - 2*rho.*((x - mx)./sx).*((y - my)./sy) ...
        + ((y - my)./sy).^2;
    p = 1./(2*pi*sx*sy*sqrt(one_minus_rho2)) ...
        .* exp(-z./(2*one_minus_rho2));
end

function gx = gaussian_pdf_1d(x, mu, sigma)
    gx = 1./(sqrt(2*pi)*sigma) .* exp(-0.5*((x - mu)./sigma).^2);
end

function F = gaussian_cdf_1d(z, mu, sigma)
    F = 0.5 .* (1 + erf((z - mu)./(sigma*sqrt(2))));
end

function MR = compute_MR(mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol)
    pdf_handle = @(x,y) gaussian_pdf_xy(x, y, mx, sx, my, sy, rho);
    MR = integral2(pdf_handle, xmin, xmax, ymin, ymax, ...
        'AbsTol', abs_tol, 'RelTol', rel_tol, 'Method', 'iterated');
end

function Pdelta = compute_Pdelta(lambda_value, delta, MR, ...
    mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    integrand = @(x) x_integrand_for_band(x, lambda_value, delta, ...
        mx, sx, my, sy, rho, ymin, ymax);

    raw_probability = integral(integrand, xmin, xmax, ...
        'ArrayValued', true, 'AbsTol', abs_tol, 'RelTol', rel_tol);

    Pdelta = raw_probability ./ MR;
    Pdelta = min(max(Pdelta, 0), 1);
end

function values = x_integrand_for_band(x, lambda_value, delta, ...
    mx, sx, my, sy, rho, ymin, ymax)

    band_lower = (x./3) .* ((1 - delta).*lambda_value - 1);
    band_upper = (x./3) .* ((1 + delta).*lambda_value - 1);

    y_lower = max(ymin, band_lower);
    y_upper = min(ymax, band_upper);

    valid_band = y_lower < y_upper;
    y_probability = zeros(size(x));

    if any(valid_band(:))
        mu_y_given_x = my + rho*(sy/sx).*(x(valid_band) - mx);
        sy_given_x = sy*sqrt(1 - rho^2);

        y_probability(valid_band) = ...
            gaussian_cdf_1d(y_upper(valid_band), mu_y_given_x, sy_given_x) ...
            - gaussian_cdf_1d(y_lower(valid_band), mu_y_given_x, sy_given_x);
    end

    values = gaussian_pdf_1d(x, mx, sx) .* y_probability;
end

function Pcurve = compute_Pcurve_for_fixed_delta(lambda_values, delta, MR, ...
    mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    Pcurve = arrayfun(@(lambda_value) compute_Pdelta(lambda_value, delta, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol), lambda_values);
end

function Qcurve = compute_Qcurve_for_fixed_lambda(delta_values, lambda_value, MR, ...
    mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

    Qcurve = arrayfun(@(delta) compute_Pdelta(lambda_value, delta, MR, ...
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol), delta_values);
end

function delta_80 = find_delta80(delta_grid, Q_grid, target_probability, ...
    lambda_peak, MR, mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol)

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
        mx, sx, my, sy, rho, xmin, xmax, ymin, ymax, abs_tol, rel_tol) ...
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
