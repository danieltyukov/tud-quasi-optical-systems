%% Q4: Decompose the oblique pattern into array factor and element pattern.
% Shows why the grating lobe at theta = -90 deg is suppressed in the
% E-plane (phi = 0) for TM incidence: the array factor has a clear lobe
% there, but the element pattern has a zero (cos^2(theta) factor in Re{Gxx}
% multiplied by the k_z = k_0 cos(theta) factor in the radiation integral
% gives an overall cos^2 dependence in voltage, cos^4 in power, that kills
% the dipole radiation along its own axis).
clear; close all; clc;
set_plot_defaults();

%% Parameters
freq    = 10e9;
c       = 3e8;
lambda  = c/freq;
k0      = 2*pi/lambda;
dx      = 20e-3;
dy      = 20e-3;
w       = 1e-3;
l       = 15e-3;
Mmax    = 30;
Nx      = 10;
Ny      = 10;
phi0    = 0;
V_TM    = 1; V_TE = 0;

theta_in = asind(lambda/dx - 1);    % 30 deg
th0      = deg2rad(theta_in);
kx0      = k0*sin(th0);
ky0      = 0;

%% Sweep
theta_obs_deg = linspace(-89.9, 89.9, 3601);
theta_obs     = deg2rad(theta_obs_deg);
kx = k0*sin(theta_obs);
ky = zeros(size(kx));
kz = k0*cos(theta_obs);

%% Element pattern: jk_z * G_xx(kx,0) * I(kx) * Jt(0)
G   = EJ_SGF(1, k0, kx, ky);
Ikx = basis_long(kx, k0, l);
Jky = basis_trans(ky, w);
EF  = (1j*kz) .* G.Gxx .* Ikx .* Jky;

%% Array factor (closed-form Dirichlet kernel)
psi_x = (kx - kx0) * dx / 2;
AFx = ones(size(psi_x)) * Nx;
nz  = abs(sin(psi_x)) > 1e-12;
AFx(nz) = sin(Nx*psi_x(nz)) ./ sin(psi_x(nz)) .* exp(1j*(Nx-1)*psi_x(nz));
psi_y = (ky - ky0) * dy / 2;
AFy = ones(size(psi_y)) * Ny;
nz  = abs(sin(psi_y)) > 1e-12;
AFy(nz) = sin(Ny*psi_y(nz)) ./ sin(psi_y(nz)) .* exp(1j*(Ny-1)*psi_y(nz));
AF = AFx .* AFy;

%% i_BF for the FSS at theta_in
Z   = Z_FSS(th0, phi0, k0, dx, dy, w, l, Mmax);
v   = v_FSS(th0, phi0, V_TM, V_TE, k0, l, w);
iBF = v / Z;

%% Total
E_tot = EF .* AF * iBF;

%% Normalize each on its own peak
P_AF  = abs(AF).^2;     P_AF_dB  = 10*log10(P_AF/max(P_AF));
P_EF  = abs(EF).^2;     P_EF_dB  = 10*log10(P_EF/max(P_EF));
P_TOT = abs(E_tot).^2;  P_TOT_dB = 10*log10(P_TOT/max(P_TOT));

%% Figure
fig = figure('Color','w','Position',[80 80 1200 460]);
set(fig, 'InvertHardcopy', 'off');
hold on; grid on;
plot(theta_obs_deg, P_AF_dB,  '-',  'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
plot(theta_obs_deg, P_EF_dB,  '--', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(theta_obs_deg, P_TOT_dB, '-',  'LineWidth', 2.0, 'Color', [0.85 0.33 0.10]);

% Theoretical GL location
s_gl  = sin(th0) - lambda/dx;
th_gl = asind(s_gl);
xline(th_gl,         'r:', 'LineWidth', 1.2);
xline(rad2deg(th0),  'k:', 'LineWidth', 1.0, 'Alpha', 0.5);
text(th_gl+1.5, -2,  sprintf('GL position\n\\theta = %+.0f°', th_gl), ...
     'FontSize', 9, 'Color', [0.7 0 0]);
text(rad2deg(th0)+1.5, -2, sprintf('main lobe\n\\theta = %+.0f°', rad2deg(th0)), ...
     'FontSize', 9, 'Color', 'k');

xlabel('\theta_{obs}  [deg]'); ylabel('Normalised |.|^2  [dB]');
title('Pattern decomposition at \theta_{inc} = 30°: array factor (grey) shows the GL at -90°, element pattern (blue) has a zero there, total (red) inherits the suppression');
legend({'Array factor |AF|^2', 'Element pattern |EF|^2', ...
        'Total pattern |EF \cdot AF|^2'}, 'Location', 'south', 'NumColumns', 3);
xlim([-90 90]); ylim([-80 1]); xticks(-90:30:90);

print('-dpng', '-r150', fullfile('..','figures','Q4_GL_suppression.png'));

fprintf('\nLevels at theoretical GL position theta = %.2f deg:\n', th_gl);
[~, iGL] = min(abs(theta_obs_deg - th_gl));
fprintf('  Array factor level (rel. to AF max):    %.2f dB\n', P_AF_dB(iGL));
fprintf('  Element pattern level (rel. to EF max): %.2f dB\n', P_EF_dB(iGL));
fprintf('  Total pattern level (rel. to TOT max):  %.2f dB\n', P_TOT_dB(iGL));
fprintf('Figure saved.\n');
