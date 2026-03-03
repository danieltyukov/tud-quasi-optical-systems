%% Q3: Reflector Far Field Patterns (2 points)

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

%% Reflector geometry for f/D = 1
fD = 1;
D_refl = f_focal / fD;
a_refl = D_refl / 2;
theta_0 = 2 * atan(D_refl / (4 * f_focal));

fprintf('f/D = %.1f, D = %.2f m, a = %.2f m\n', fD, D_refl, a_refl);
fprintf('theta_0 = %.2f deg\n', theta_0 * 180/pi);
fprintf('First Airy null at: %.5f deg\n', asin(1.22*lambda0/D_refl)*180/pi);

%% Compute aperture currents on grid
N = 301;
x_ap = linspace(-a_refl, a_refl, N);
y_ap = linspace(-a_refl, a_refl, N);
dx = x_ap(2) - x_ap(1);
dy = y_ap(2) - y_ap(1);
[X_ap, Y_ap] = meshgrid(x_ap, y_ap);
rho_ap = sqrt(X_ap.^2 + Y_ap.^2);
mask = rho_ap <= a_refl;

theta_p = 2 * atan(rho_ap / (2 * f_focal));
phi_p = atan2(Y_ap, X_ap);
phi_p(rho_ap < 1e-15) = 0;

[E_a_x, E_a_y] = aperture_field(k0, f_focal, a_feed, theta_p, phi_p);
E_a_x = E_a_x .* mask;
E_a_y = E_a_y .* mask;

M_x = E_a_y;
M_y = -E_a_x;
J_x = -E_a_x / zeta;
J_y = -E_a_y / zeta;

%% Far field via numerical FT (Balanis spatial domain)
theta_null_airy = asin(1.22 * lambda0 / D_refl);
theta_max = 15 * theta_null_airy;
Nth_ff = 2001;
theta_ff = linspace(0, theta_max, Nth_ff);

planes = struct('phi', {pi/2, 0}, 'label', {'E-plane (\phi=90°)', 'H-plane (\phi=0°)'});
E_tot_E = zeros(1, Nth_ff);
E_tot_H = zeros(1, Nth_ff);

for ip = 1:2
    phi_ff = planes(ip).phi;
    for it = 1:Nth_ff
        th = theta_ff(it);
        kx_ff = k0 * sin(th) * cos(phi_ff);
        ky_ff = k0 * sin(th) * sin(phi_ff);

        phase = exp(1j * (kx_ff * X_ap + ky_ff * Y_ap));
        Nx = sum(sum(J_x .* phase)) * dx * dy;
        Ny = sum(sum(J_y .* phase)) * dx * dy;
        Lx = sum(sum(M_x .* phase)) * dx * dy;
        Ly = sum(sum(M_y .* phase)) * dx * dy;

        cos_th = cos(th);
        cos_ph = cos(phi_ff);
        sin_ph = sin(phi_ff);

        % Balanis: E_th ~ L_phi + zeta*N_theta, E_ph ~ L_theta - zeta*N_phi
        E_th = (-Lx * sin_ph + Ly * cos_ph) + zeta * (Nx * cos_ph + Ny * sin_ph) * cos_th;
        E_ph = (Lx * cos_ph + Ly * sin_ph) * cos_th + zeta * (Nx * sin_ph - Ny * cos_ph);

        if ip == 1
            E_tot_E(it) = sqrt(abs(E_th)^2 + abs(E_ph)^2);
        else
            E_tot_H(it) = sqrt(abs(E_th)^2 + abs(E_ph)^2);
        end
    end
end

E_max = max([max(E_tot_E), max(E_tot_H)]);
E_E_dB = 20 * log10(E_tot_E / E_max);
E_H_dB = 20 * log10(E_tot_H / E_max);

%% Airy pattern reference
airy = ones(size(theta_ff));
for it = 1:Nth_ff
    arg = k0 * a_refl * sin(theta_ff(it));
    if abs(arg) > 1e-10
        airy(it) = 2 * besselj(1, arg) / arg;
    end
end
airy_dB = 20 * log10(abs(airy) / max(abs(airy)));

%% Plot
fig = figure('Color', 'w', 'Position', [100 100 900 550]);
set(fig, 'InvertHardcopy', 'off');
plot(theta_ff * 180/pi, E_E_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_ff * 180/pi, E_H_dB, 'r--', 'LineWidth', 1.5);
plot(theta_ff * 180/pi, airy_dB, 'k:', 'LineWidth', 1.5);
xlabel('\theta [deg]');
ylabel('Normalized |E^{far}| [dB]');
title('Q3: Reflector Far Field (f/D = 1, D = 5 m, 500 GHz)');
legend('E-plane (\phi = 90°)', 'H-plane (\phi = 0°)', 'Airy pattern (uniform)', ...
       'Location', 'northeast');
grid on;
ylim([-40 0]);
xlim([0 theta_max * 180/pi]);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q3_reflector_farfield.png'));
