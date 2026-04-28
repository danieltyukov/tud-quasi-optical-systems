%% Q2: scan analysis of the chosen array.
%% Compute |Gamma_act(theta_scan, f)| in E and H planes; locate GL onset
%% and any scan-blindness behaviour; compute the maximum usable scan angle.

clear; close all; clc;
set_plot_defaults();
load('Q1_design.mat', 'h_best','l_best','w','dx','dy','ZL_best','Mmax', ...
                      'f_lo','f_hi','f_c','lam_c');

c0 = 3e8;

% scan-angle grid
theta_deg = linspace(-89, 89, 359);
theta     = deg2rad(theta_deg);

% Frequency grid: band edges + center
freqs = [f_lo, f_c, f_hi];

% Three principal scan planes
planes = struct( ...
    'name',  {'E-plane (\phi_0 = 0\circ)', 'H-plane (\phi_0 = 90\circ)', 'D-plane (\phi_0 = 45\circ)'}, ...
    'short', {'E','H','D'}, ...
    'phi',   {0, pi/2, pi/4});

GammaDB_grid = nan(numel(freqs), numel(planes), numel(theta));
Z_grid       = complex(nan(numel(freqs), numel(planes), numel(theta)));

fprintf('Computing |Gamma_act(theta, f)|...\n');
for kf = 1:numel(freqs)
    k0 = 2*pi*freqs(kf)/c0;
    for kp = 1:numel(planes)
        for kt = 1:numel(theta)
            Z = Z_active_GP(theta(kt), planes(kp).phi, k0, dx, dy, w, l_best, h_best, Mmax);
            Z_grid(kf, kp, kt) = Z;
            GammaDB_grid(kf, kp, kt) = 20*log10(abs((Z - ZL_best)/(Z + ZL_best)));
        end
    end
    fprintf('  f = %.2f GHz done (%d/%d)\n', freqs(kf)/1e9, kf, numel(freqs));
end

% Theoretical grating-lobe entry angle for each frequency
sin_GL = c0./freqs/dx - 1;             % might be > 1 (no GL anywhere)
theta_GL_deg = nan(size(freqs));
for kf = 1:numel(freqs)
    if sin_GL(kf) <= 1
        theta_GL_deg(kf) = asind(sin_GL(kf));
    end
end

% Find max scan angle satisfying |Gamma| < -10 dB at all three frequencies
% in each plane (worst-case across the band)
worst_band = squeeze(max(GammaDB_grid, [], 1));     % [planes x theta]
% Allow up to where worst_band crosses -10 dB
theta_max = nan(numel(planes),1);
for kp = 1:numel(planes)
    g = squeeze(worst_band(kp, :));
    % Take the contiguous region around theta=0 where g < -10 dB
    ok = g < -10;
    [~, i0] = min(abs(theta_deg));
    if ~ok(i0)
        theta_max(kp) = 0;
        continue;
    end
    iL = i0; while iL > 1 && ok(iL-1), iL = iL-1; end
    iR = i0; while iR < numel(theta_deg) && ok(iR+1), iR = iR+1; end
    theta_max(kp) = min(theta_deg(iR), -theta_deg(iL));
end

fprintf('\n=== Maximum scan angle satisfying |Gamma| < -10 dB across band ===\n');
for kp = 1:numel(planes)
    fprintf('  %s : theta_max = %5.1f deg\n', planes(kp).short, theta_max(kp));
end
for kf = 1:numel(freqs)
    if isnan(theta_GL_deg(kf))
        fprintf('  GL onset at f = %.2f GHz : no GL in visible region (lambda/d - 1 = %.3f)\n', ...
                freqs(kf)/1e9, sin_GL(kf));
    else
        fprintf('  GL onset at f = %.2f GHz : theta_GL = %.2f deg\n', ...
                freqs(kf)/1e9, theta_GL_deg(kf));
    end
end

%% Plot 1 : |Gamma_act| vs scan angle, per plane, colored by frequency
fig = figure('Color','w','Position',[60 60 1300 400]);
set(fig, 'InvertHardcopy', 'off');

f_colors  = [0.0 0.4 0.7; 0.0 0.0 0.0; 0.85 0.33 0.10];
f_labels  = {sprintf('27.5 GHz'), sprintf('29.25 GHz'), sprintf('31.0 GHz')};

for kp = 1:numel(planes)
    subplot(1, numel(planes), kp); hold on; grid on;
    for kf = 1:numel(freqs)
        plot(theta_deg, squeeze(GammaDB_grid(kf, kp, :)), '-', ...
             'Color', f_colors(kf,:), 'LineWidth', 1.5);
    end
    yline(-10, 'k:', 'LineWidth', 1);
    xline(theta_max(kp), 'r--', sprintf('%.0f\\circ', theta_max(kp)), ...
          'LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
    xline(-theta_max(kp), 'r--', sprintf('-%.0f\\circ', theta_max(kp)), ...
          'LabelVerticalAlignment','top','LabelHorizontalAlignment','right');
    if kp==1, legend(f_labels, 'Location', 'south'); end
    xlabel('\theta_0 [deg]'); ylabel('|\Gamma_{act}| [dB]');
    xlim([-90 90]); ylim([-30 0]); xticks(-90:30:90);
    title(planes(kp).name);
end
sgtitle({'Q2: active reflection coefficient vs scan angle, three frequencies in the band', ...
         sprintf('h = %.2f mm, l = %.2f mm, w = %.2f mm, d_x = d_y = %.2f mm, Z_L = %d \\Omega', ...
                 h_best*1e3, l_best*1e3, w*1e3, dx*1e3, ZL_best)}, ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng','-r150', fullfile('..','figures','Q2_Gamma_vs_scan.png'));

%% Plot 2 : 2D map of |Gamma| in (f, theta) plane, H-plane (worst case)
fine_freqs = linspace(0.92*f_lo, 1.05*f_hi, 60);
fine_theta_deg = linspace(-89, 89, 181);
GammaH_2D = nan(numel(fine_freqs), numel(fine_theta_deg));
fprintf('\nFine 2D scan map (H-plane)...\n');
for kf = 1:numel(fine_freqs)
    k0 = 2*pi*fine_freqs(kf)/c0;
    for kt = 1:numel(fine_theta_deg)
        th = deg2rad(fine_theta_deg(kt));
        Z  = Z_active_GP(th, pi/2, k0, dx, dy, w, l_best, h_best, Mmax);
        GammaH_2D(kf, kt) = 20*log10(abs((Z - ZL_best)/(Z + ZL_best)));
    end
end

fig2 = figure('Color','w','Position',[60 60 850 500]);
set(fig2, 'InvertHardcopy', 'off');
GammaH_2D_clipped = max(min(GammaH_2D, 0), -30);
imagesc(fine_theta_deg, fine_freqs/1e9, GammaH_2D_clipped);
set(gca, 'YDir', 'normal');
colormap(parula);
cb = colorbar; ylabel(cb, '|\Gamma_{act}| [dB]');
caxis([-30 0]);
hold on;
contour(fine_theta_deg, fine_freqs/1e9, GammaH_2D_clipped, [-10 -10], 'r-', 'LineWidth', 2);
yline(f_lo/1e9, 'w-', 'LineWidth', 1.4); yline(f_hi/1e9, 'w-', 'LineWidth', 1.4);
xlabel('\theta_0 [deg]'); ylabel('frequency [GHz]');
title('H-plane |\Gamma_{act}|: red contour is the -10 dB usable region');
print('-dpng','-r150', fullfile('..','figures','Q2_Gamma_2D_Hplane.png'));

%% Plot 3 : grating-lobe diagram for the chosen lattice at f_lo, f_c, f_hi
fig3 = figure('Color','w','Position',[60 60 1100 380]);
set(fig3, 'InvertHardcopy', 'off');
for kf = 1:numel(freqs)
    k0    = 2*pi*freqs(kf)/c0;
    lambda = c0/freqs(kf);
    rho_x = lambda/dx;
    rho_y = lambda/dy;
    subplot(1,3,kf); hold on; axis equal;
    th_circ = linspace(0,2*pi,200);
    % Visible region (radius 1 in normalized k-space)
    fill(cos(th_circ), sin(th_circ), [0.85 0.92 1.0], 'EdgeColor', [0.2 0.4 0.7], 'LineWidth', 1.4);
    % Higher-order Floquet circles (centred at integer multiples of (lambda/d))
    for mx = -2:2
        for my = -2:2
            if mx==0 && my==0, continue; end
            cx = mx*rho_x; cy = my*rho_y;
            if (abs(cx) <= 2.5) && (abs(cy) <= 2.5)
                plot(cx + cos(th_circ), cy + sin(th_circ), '-', 'Color', [0.5 0.5 0.5]);
                text(cx, cy, sprintf('(%d,%d)', mx, my), 'HorizontalAlignment','center', ...
                     'FontSize', 7, 'Color', [0.3 0.3 0.3]);
            end
        end
    end
    % E-plane and H-plane scan trajectories
    plot([-1 1], [0 0], 'r--', 'LineWidth', 1.2);
    plot([0 0], [-1 1], 'g--', 'LineWidth', 1.2);
    grid on;
    xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    xlabel('k_{x0}/k_0'); ylabel('k_{y0}/k_0');
    title(sprintf('f = %.2f GHz   (\\lambda/d = %.3f)', freqs(kf)/1e9, lambda/dx));
end
sgtitle('Q2: grating-lobe diagram for the chosen lattice (d = 4.8 mm)', 'FontSize', 11, 'FontWeight', 'bold');
print('-dpng','-r150', fullfile('..','figures','Q2_GLdiagram.png'));

save('Q2_scan.mat', 'theta_deg','freqs','planes','GammaDB_grid','Z_grid', ...
     'theta_GL_deg','theta_max','sin_GL','fine_freqs','fine_theta_deg','GammaH_2D');
fprintf('\nFigures saved.\n');
