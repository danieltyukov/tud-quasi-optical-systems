%% Q2: Aperture Current Distribution (2 points)

clear; close all; clc;
set_plot_defaults();

%% Constants
c = 3e8;
freq = 500e9;
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
zeta = 120 * pi;
D_f = 4 * lambda0;
a_feed = D_f / 2;
f_focal = 5;

%% Reflector geometry for f/D = 2
fD = 2;
D_refl = f_focal / fD;
a_refl = D_refl / 2;
theta_0 = 2 * atan(D_refl / (4 * f_focal));
fprintf('f/D = %.1f, D = %.2f m, theta_0 = %.2f deg\n', fD, D_refl, theta_0*180/pi);

%% Aperture grid
N = 201;
x = linspace(-a_refl, a_refl, N);
y = linspace(-a_refl, a_refl, N);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);
mask = rho <= a_refl;

theta_p = 2 * atan(rho / (2 * f_focal));
phi_p = atan2(Y, X);
phi_p(rho < 1e-15) = 0;

%% Compute aperture field and equivalent currents
[E_a_x, E_a_y] = aperture_field(k0, f_focal, a_feed, theta_p, phi_p);
E_a_x = E_a_x .* mask;
E_a_y = E_a_y .* mask;

M_x = E_a_y;           % co-pol
M_y = -E_a_x;          % cross-pol
J_x = -E_a_x / zeta;   % cross-pol
J_y = -E_a_y / zeta;    % co-pol

%% x-polarization figure
fig1 = figure('Color', 'w', 'Position', [50 50 900 400]);
set(fig1, 'InvertHardcopy', 'off');

ax1 = subplot(1,2,1);
Mx_dB = 20*log10(abs(M_x) / max(abs(M_x(mask))));
Mx_dB(~mask) = NaN;
imagesc(x*1e3, y*1e3, Mx_dB);
axis equal tight; set(gca, 'YDir', 'normal', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
colormap(ax1, 'jet'); clim([-40 0]); colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('|M_{s,x}| [dB] (co-pol)');

ax2 = subplot(1,2,2);
Jx_dB = 20*log10(abs(J_x) / max(abs(J_x(mask))));
Jx_dB(~mask) = NaN;
imagesc(x*1e3, y*1e3, Jx_dB);
axis equal tight; set(gca, 'YDir', 'normal', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
colormap(ax2, 'jet'); clim([-40 0]); colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('|J_{s,x}| [dB] (cross-pol)');

sgtitle('Q2: x-polarization Aperture Currents (f/D = 2)', 'Color', 'k', 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_aperture_currents_xpol.png'));

%% y-polarization figure
fig2 = figure('Color', 'w', 'Position', [50 50 900 400]);
set(fig2, 'InvertHardcopy', 'off');

ax3 = subplot(1,2,1);
My_dB = 20*log10(abs(M_y) / max(abs(M_y(mask))));
My_dB(~mask) = NaN;
imagesc(x*1e3, y*1e3, My_dB);
axis equal tight; set(gca, 'YDir', 'normal', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
colormap(ax3, 'jet'); clim([-40 0]); colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('|M_{s,y}| [dB] (cross-pol)');

ax4 = subplot(1,2,2);
Jy_dB = 20*log10(abs(J_y) / max(abs(J_y(mask))));
Jy_dB(~mask) = NaN;
imagesc(x*1e3, y*1e3, Jy_dB);
axis equal tight; set(gca, 'YDir', 'normal', 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
colormap(ax4, 'jet'); clim([-40 0]); colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('|J_{s,y}| [dB] (co-pol)');

sgtitle('Q2: y-polarization Aperture Currents (f/D = 2)', 'Color', 'k', 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_aperture_currents_ypol.png'));
