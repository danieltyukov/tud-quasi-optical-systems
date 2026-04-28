%% Q1: Active input impedance vs scan angle for E, H, D planes.
% Two unit-cell sizes (15 mm and 20 mm), printed sinusoidal-current dipoles
% (l = 14 mm, w = 1 mm) on a free-space rectangular grid at 10 GHz.
clear; close all; clc;
set_plot_defaults();

%% Parameters
freq    = 10e9;
c       = 3e8;
lambda  = c/freq;
k0      = 2*pi/lambda;
w       = 1e-3;
l       = 14e-3;
Mmax    = 30;          % Floquet truncation, |m| <= Mmax in each direction

cases = struct( ...
    'name',  {'15 mm', '20 mm'}, ...
    'tag',   {'15mm',  '20mm'}, ...
    'dx',    {15e-3,   20e-3}, ...
    'dy',    {15e-3,   20e-3});

% Scan angles: avoid exact endfire to keep the fundamental kz=0 spike finite
theta_deg = linspace(-89, 89, 1441);
theta     = deg2rad(theta_deg);

planes = struct( ...
    'name',  {'E-plane (\phi = 0\circ)', 'H-plane (\phi = 90\circ)', 'D-plane (\phi = 45\circ)'}, ...
    'short', {'E', 'H', 'D'}, ...
    'phi',   {0, pi/2, pi/4});

plane_styles = {'-', '--', '-.'};
plane_colors = [0 0.45 0.74;        % E   blue
                0.85 0.33 0.10;     % H   orange-red
                0.47 0.67 0.19];    % D   green

%% Compute Z_in for every (case, plane, theta)
fprintf('Computing Z_in...\n');
results = cell(numel(cases), numel(planes));
for ci = 1:numel(cases)
    for pi_ = 1:numel(planes)
        Z = zeros(size(theta));
        for ti = 1:numel(theta)
            Z(ti) = Z_active(theta(ti), planes(pi_).phi, k0, ...
                             cases(ci).dx, cases(ci).dy, w, l, Mmax);
        end
        results{ci, pi_} = Z;
        fprintf('  %s, %s done\n', cases(ci).name, planes(pi_).short);
    end
end

%% --- Figure 1: Re/Im of Z_in vs theta, two rows (cases) x two cols (Re, Im) ---
fig = figure('Color','w','Position',[80 80 1200 850]);
set(fig, 'InvertHardcopy', 'off');

for ci = 1:numel(cases)
    % --- Resistance ---
    subplot(numel(cases), 2, 2*(ci-1)+1); hold on; grid on;
    h = gobjects(1, numel(planes));
    for pi_ = 1:numel(planes)
        h(pi_) = plot(theta_deg, real(results{ci,pi_}), plane_styles{pi_}, ...
                      'LineWidth', 1.6, 'Color', plane_colors(pi_,:));
    end
    xlabel('\theta_0 [deg]'); ylabel('Re\{Z_{in}\} [\Omega]');
    title(sprintf('d_x = d_y = %s   —   Re\\{Z_{in}\\}', cases(ci).name));
    xlim([-90 90]); ylim([0 350]); xticks(-90:30:90);
    if ci == 1
        legend(h, {planes.name}, 'Location', 'north', 'NumColumns', 3);
    end

    % --- Reactance ---
    subplot(numel(cases), 2, 2*(ci-1)+2); hold on; grid on;
    for pi_ = 1:numel(planes)
        plot(theta_deg, imag(results{ci,pi_}), plane_styles{pi_}, ...
             'LineWidth', 1.6, 'Color', plane_colors(pi_,:));
    end
    xlabel('\theta_0 [deg]'); ylabel('Im\{Z_{in}\} [\Omega]');
    title(sprintf('d_x = d_y = %s   —   Im\\{Z_{in}\\}', cases(ci).name));
    xlim([-90 90]); ylim([-350 350]); xticks(-90:30:90);

    % Mark the H-plane GL entry for the 20 mm case
    if ci == 2
        % theta_GL = asind(lambda/dy - 1)
        thG = asind(lambda/cases(ci).dy - 1);
        subplot(numel(cases), 2, 2*(ci-1)+1);
        xline(+thG, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
        xline(-thG, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
        subplot(numel(cases), 2, 2*(ci-1)+2);
        xline(+thG, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
        xline(-thG, 'k:', 'LineWidth', 1, 'Alpha', 0.5);
    end
end

sgtitle({'Active input impedance of an infinite array of printed dipoles', ...
         sprintf('w = %g mm,   l = %g mm,   f = %g GHz   (\\lambda_0 = %g mm)', ...
                 w*1e3, l*1e3, freq/1e9, lambda*1e3)}, ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q1_Zin_vs_angle.png'));

%% --- Figure 2: |Z_in| on log scale -- captures all singularities at once ---
fig2 = figure('Color','w','Position',[80 80 1200 450]);
set(fig2, 'InvertHardcopy', 'off');
for ci = 1:numel(cases)
    subplot(1, 2, ci); hold on; grid on;
    for pi_ = 1:numel(planes)
        semilogy(theta_deg, abs(results{ci,pi_}), plane_styles{pi_}, ...
                 'LineWidth', 1.6, 'Color', plane_colors(pi_,:));
    end
    set(gca, 'YScale', 'log');
    xlabel('\theta_0 [deg]'); ylabel('|Z_{in}| [\Omega]');
    title(sprintf('d_x = d_y = %s', cases(ci).name));
    xlim([-90 90]); ylim([10 1e4]); xticks(-90:30:90);
    if ci == 2
        thG = asind(lambda/cases(ci).dy - 1);
        xline(+thG, 'k:', 'LineWidth', 1);
        xline(-thG, 'k:', 'LineWidth', 1);
        text(thG+1, 5e3, '\theta_{GL} = 30°', 'FontSize', 9);
    end
    if ci == 1
        legend({planes.name}, 'Location', 'north', 'NumColumns', 3);
    end
end
sgtitle('|Z_{in}| vs scan angle (log scale): singularities at the H-plane GL entry and at endfire', ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q1_Zin_magnitude.png'));

%% --- Figure 3: zoom on the H-plane scan-blindness for 20 mm ---
fig3 = figure('Color','w','Position',[80 80 1100 400]);
set(fig3, 'InvertHardcopy', 'off');
ZH20 = results{2, 2};
thG = asind(lambda/cases(2).dy - 1);

subplot(1,2,1); hold on; grid on;
plot(theta_deg, real(ZH20), 'Color', plane_colors(2,:), 'LineWidth', 1.7);
xline(+thG, 'k:', 'LineWidth', 1); xline(-thG, 'k:', 'LineWidth', 1);
xlabel('\theta_0 [deg]'); ylabel('Re\{Z_{in}\} [\Omega]');
title('20 mm — H-plane resistance');
xlim([-90 90]); ylim([0 700]); xticks(-90:30:90);
text(31, 640, '\theta_0 = +30°', 'FontSize', 9);

subplot(1,2,2); hold on; grid on;
plot(theta_deg, imag(ZH20), 'Color', plane_colors(2,:), 'LineWidth', 1.7);
xline(+thG, 'k:', 'LineWidth', 1); xline(-thG, 'k:', 'LineWidth', 1);
xlabel('\theta_0 [deg]'); ylabel('Im\{Z_{in}\} [\Omega]');
title('20 mm — H-plane reactance');
xlim([-90 90]); ylim([-700 700]); xticks(-90:30:90);

sgtitle('Scan blindness in the H-plane: Floquet mode (m_y = \pm 1) enters at \theta_0 = \pm 30°', ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q1_Zin_zoom_Hplane_20mm.png'));

%% --- Save numeric data for the report ---
broadside = struct();
for ci = 1:numel(cases)
    Z0 = Z_active(0, 0, k0, cases(ci).dx, cases(ci).dy, w, l, Mmax);
    broadside.(['c' cases(ci).tag]) = Z0;
end
save('Q1_results.mat', 'theta_deg','results','cases','planes','broadside','freq','k0','lambda');

fprintf('\nBroadside Z_in:\n');
fprintf('  15 mm: Z = %.3f + j(%.3f) Ohm\n', real(broadside.c15mm), imag(broadside.c15mm));
fprintf('  20 mm: Z = %.3f + j(%.3f) Ohm\n', real(broadside.c20mm), imag(broadside.c20mm));

% Numeric peak in Re for 20 mm H-plane (clipped, for comment in report)
[mx_, ix] = max(real(results{2,2}));
fprintf('\n20 mm H-plane peak Re{Z}: %.1f Ohm at theta = %.2f deg (numerical, grid-limited)\n', ...
        mx_, theta_deg(ix));
fprintf('Theoretical GL entry angle: %.4f deg\n', asind(lambda/cases(2).dy - 1));
fprintf('\nFigures saved.\n');
