%% Q1: Active element pattern vs incidence angle, TM, phi = 0 plane.
clear; close all; clc;
set_plot_defaults();

freq    = 10e9;
c       = 3e8;
lambda  = c/freq;
k0      = 2*pi/lambda;
dx      = 20e-3;
dy      = 20e-3;
w       = 1e-3;
l       = 15e-3;
Mmax    = 30;
phi0    = 0;
V_TM    = 1; V_TE = 0;

sin_in   = lambda/dx - 1;
theta_in = asind(sin_in);
fprintf('Theoretical GL entry: sin(theta_in) = lambda/d - 1 = %.4f\n', sin_in);
fprintf('                      theta_in       = %.4f deg\n', theta_in);

theta_deg = linspace(-89, 89, 1791);
theta     = deg2rad(theta_deg);

Z   = zeros(size(theta));
v   = zeros(size(theta));
iBF = zeros(size(theta));
EP  = zeros(size(theta));
AEP = zeros(size(theta));

for ti = 1:numel(theta)
    th = theta(ti);
    kx0 = k0*sin(th)*cos(phi0);
    ky0 = k0*sin(th)*sin(phi0);

    Z(ti)   = Z_FSS(th, phi0, k0, dx, dy, w, l, Mmax);
    v(ti)   = v_FSS(th, phi0, V_TM, V_TE, k0, l, w);
    iBF(ti) = v(ti) / Z(ti);

    G    = EJ_SGF(1, k0, kx0, ky0);
    Ikx0 = basis_long(kx0, k0, l);
    Jky0 = basis_trans(ky0, w);
    kz0  = k0*cos(th);

    EP(ti)  = 1j*kz0 * G.Gxx * Ikx0 * Jky0;
    AEP(ti) = EP(ti) * iBF(ti);
end

[~, i0]  = min(abs(theta_deg));
fprintf('\nBroadside (theta_0 = 0):\n');
fprintf('  Z        = %.3f + j(%.3f) Ohm\n', real(Z(i0)), imag(Z(i0)));
fprintf('  v        = %.4e\n',                 abs(v(i0)));
fprintf('  i_BF     = %.4e\n',                 abs(iBF(i0)));
fprintf('  |AEP|^2  = %.4e (linear, peak ref)\n', abs(AEP(i0))^2);

P_AEP    = abs(AEP).^2;
P_AEP_dB = 10*log10(P_AEP / max(P_AEP));

fig = figure('Color','w','Position',[80 80 1100 460]);
set(fig, 'InvertHardcopy', 'off');

subplot(1,2,1); hold on; grid on;
plot(theta_deg, P_AEP_dB, 'LineWidth', 1.7, 'Color', [0 0.45 0.74]);
xline(+theta_in, 'r:', 'LineWidth', 1.2);
xline(-theta_in, 'r:', 'LineWidth', 1.2);
text(+theta_in+1.0, -1.5, sprintf('\\theta_{in} = %.0f°', theta_in), ...
     'FontSize', 9, 'Color', [0.7 0 0]);
xlabel('\theta_0  [deg]'); ylabel('|AEP|^2 / |AEP|^2_{max}  [dB]');
title('Active element pattern (dB)');
xlim([-90 90]); ylim([-30 1]); xticks(-90:30:90);

subplot(1,2,2); hold on; grid on;
plot(theta_deg, real(Z), '-',  'LineWidth', 1.6, 'Color', [0 0.45 0.74]);
plot(theta_deg, imag(Z), '--', 'LineWidth', 1.6, 'Color', [0.85 0.33 0.10]);
xline(+theta_in, 'r:', 'LineWidth', 1.2);
xline(-theta_in, 'r:', 'LineWidth', 1.2);
xlabel('\theta_0  [deg]'); ylabel('Z_{FSS}  [\Omega]');
title('FSS self-impedance Z = -\Sigma G_{xx} |B|^2 / (d_x d_y)');
legend('Re\{Z\}', 'Im\{Z\}', 'Location', 'best');
xlim([-90 90]); ylim([-150 350]); xticks(-90:30:90);

sgtitle({'FSS active element pattern -- TM incidence, \phi_0 = 0 plane', ...
         sprintf('d_x = d_y = %g mm,   l = %g mm,   w = %g mm,   f = %g GHz', ...
                 dx*1e3, l*1e3, w*1e3, freq/1e9)}, ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q1_AEP_vs_theta.png'));

save('Q1_AEP_results.mat', 'theta_deg','Z','v','iBF','EP','AEP', ...
     'theta_in','sin_in','freq','k0','lambda','dx','dy','w','l');

fprintf('\nFigure saved.\n');
