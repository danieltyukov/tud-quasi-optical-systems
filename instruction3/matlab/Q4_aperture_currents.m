%% Q4: Aperture Current Distribution

clear; close all; clc;
set_plot_defaults();

c = 3e8; freq = 200e9; lambda0 = c / freq;
k0 = 2*pi/lambda0; er = 4; zeta_0 = 120*pi;
D = 9*lambda0; theta_0 = 40*pi/180;
u0 = 0.5; v0 = 0.5;

lens = ellipse_geometry(er, D, lambda0);
e = lens.e; a = lens.a;

Npts = 501;
theta_arr = linspace(0, theta_0, Npts);
r_ell = a*(1-e^2) ./ (1 - e*cos(theta_arr));
rho = r_ell .* sin(theta_arr);
theta_i = lens_incidence_angle(e, theta_arr);
[tau_perp, tau_par, ~, ~] = fresnel_coeff(er, theta_i);

% Helper to compute Jx, Jy for a given phi
compute_J = @(phi_val) deal( ...
    -(2/zeta_0) * (tau_par .* (sin(phi_val)*exp(-((sin(theta_arr)*cos(phi_val)/u0).^2 + (sin(theta_arr)*sin(phi_val)/v0).^2))) .* cos(phi_val) ...
    - tau_perp .* (cos(phi_val)*exp(-((sin(theta_arr)*cos(phi_val)/u0).^2 + (sin(theta_arr)*sin(phi_val)/v0).^2))) .* sin(phi_val)) ./ r_ell, ...
    -(2/zeta_0) * (tau_par .* (sin(phi_val)*exp(-((sin(theta_arr)*cos(phi_val)/u0).^2 + (sin(theta_arr)*sin(phi_val)/v0).^2))) .* sin(phi_val) ...
    + tau_perp .* (cos(phi_val)*exp(-((sin(theta_arr)*cos(phi_val)/u0).^2 + (sin(theta_arr)*sin(phi_val)/v0).^2))) .* cos(phi_val)) ./ r_ell);

% E-plane (phi=90): E_theta = gauss, E_phi = 0 -> Jy uses tau_par
gauss_E = exp(-(sin(theta_arr)/v0).^2);
Jy_E = -(2/zeta_0) * tau_par .* gauss_E ./ r_ell;
Jx_E = zeros(size(theta_arr));

% H-plane (phi=0): E_theta = 0, E_phi = gauss -> Jy uses tau_perp
gauss_H = exp(-(sin(theta_arr)/u0).^2);
Jy_H = -(2/zeta_0) * tau_perp .* gauss_H ./ r_ell;
Jx_H = zeros(size(theta_arr));

% phi=45 plane: both E_theta and E_phi nonzero -> cross-pol appears
phi45 = pi/4;
gauss_45 = exp(-((sin(theta_arr)*cos(phi45)/u0).^2 + (sin(theta_arr)*sin(phi45)/v0).^2));
Ei_th_45 = sin(phi45) * gauss_45;
Ei_ph_45 = cos(phi45) * gauss_45;
Jx_45 = -(2/zeta_0) * (tau_par.*Ei_th_45*cos(phi45) - tau_perp.*Ei_ph_45*sin(phi45)) ./ r_ell;
Jy_45 = -(2/zeta_0) * (tau_par.*Ei_th_45*sin(phi45) + tau_perp.*Ei_ph_45*cos(phi45)) ./ r_ell;

J_max = max(abs(Jy_E));
rho_norm = rho / (D/2);

fprintf('Edge taper Jy E-plane (tau_par): %.1f dB\n', 20*log10(abs(Jy_E(end))/abs(Jy_E(1))));
fprintf('Edge taper Jy H-plane (tau_perp): %.1f dB\n', 20*log10(abs(Jy_H(end))/abs(Jy_H(1))));
fprintf('Max cross-pol Jx at phi=45: %.1f dB below co-pol\n', ...
    20*log10(max(abs(Jx_45))/J_max));

fig = figure('Color','w','Position',[100 100 900 650]);
set(fig,'InvertHardcopy','off');

subplot(2,1,1);
plot(rho_norm, 20*log10(abs(Jy_E)/J_max+1e-15), 'b-', 'LineWidth', 2); hold on;
plot(rho_norm, 20*log10(abs(Jy_H)/J_max+1e-15), 'r--', 'LineWidth', 2);
plot(rho_norm, 20*log10(abs(Jy_45)/J_max+1e-15), 'g-.', 'LineWidth', 1.5);
xlabel('\rho / (D/2)'); ylabel('Normalized |J_y| [dB]');
title('Q4: Co-pol Aperture Current (J_y)');
lg = legend('E-plane \phi=90° (\tau^\parallel)', 'H-plane \phi=0° (\tau^\perp)', ...
       '\phi=45°', 'Location', 'southwest');
set(lg, 'FontSize', 10);
grid on; ylim([-25 2]); xlim([0 1]);

subplot(2,1,2);
plot(rho_norm, 20*log10(abs(Jx_45)/J_max+1e-15), 'm-', 'LineWidth', 2); hold on;
plot(rho_norm, 20*log10(abs(Jx_E)/J_max+1e-15), 'b:', 'LineWidth', 1.5);
plot(rho_norm, 20*log10(abs(Jx_H)/J_max+1e-15), 'r:', 'LineWidth', 1.5);
xlabel('\rho / (D/2)'); ylabel('Normalized |J_x| [dB]');
title('Q4: Cross-pol Aperture Current (J_x)');
lg = legend('\phi=45° (cross-pol)', 'E-plane \phi=90° (zero)', ...
       'H-plane \phi=0° (zero)', 'Location', 'southwest');
set(lg, 'FontSize', 10);
grid on; ylim([-50 0]); xlim([0 1]);

print('-dpng', '-r150', fullfile('..', 'figures', 'Q4_aperture_currents.png'));
fprintf('Q4 done.\n');
