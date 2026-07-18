%% Q5: Directivity and Gain vs Diameter (2 points)

clear; close all; clc;
set_plot_defaults();

%% Constants
c = 3e8;
freq = 500e9;
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
D_f = 4 * lambda0;
a_feed = D_f / 2;
f_focal = 5;

%% Pre-compute feed pattern (same as Q4)
Nth_feed = 901;
Nph_feed = 361;
theta_feed = linspace(0, pi/2, Nth_feed);
phi_feed = linspace(0, 2*pi, Nph_feed);
dth = theta_feed(2) - theta_feed(1);
dph = phi_feed(2) - phi_feed(1);
[TH_feed, PH_feed] = meshgrid(theta_feed, phi_feed);

[E_th_2d, E_ph_2d] = feed_farfield(k0, a_feed, TH_feed, PH_feed);
U_feed = abs(E_th_2d).^2 + abs(E_ph_2d).^2;
E_co = E_th_2d .* sin(PH_feed) + E_ph_2d .* cos(PH_feed);
G_feed = U_feed;
P_total = sum(sum(U_feed .* sin(TH_feed) * dth * dph));

%% Sweep f/D
N_fD = 200;
fD_range = linspace(0.6, 6, N_fD);
D_max_dBi = zeros(1, N_fD);
D_dBi = zeros(1, N_fD);
G_dBi = zeros(1, N_fD);

for i_fD = 1:N_fD
    fD = fD_range(i_fD);
    D_refl = f_focal / fD;
    theta_0 = 2 * atan(D_refl / (4 * f_focal));
    D_max = (pi * D_refl / lambda0)^2;
    D_max_dBi(i_fD) = 10 * log10(D_max);

    idx_th0 = find(theta_feed <= theta_0, 1, 'last');
    P_intercepted = sum(sum(U_feed(:, 1:idx_th0) .* sin(TH_feed(:, 1:idx_th0)) * dth * dph));
    eta_s = P_intercepted / P_total;

    integrand_num = abs(E_co(:, 1:idx_th0)) .* tan(TH_feed(:, 1:idx_th0)/2);
    I_num = sum(sum(integrand_num * dth * dph));
    integrand_den = G_feed(:, 1:idx_th0) .* sin(TH_feed(:, 1:idx_th0));
    I_den = sum(sum(integrand_den * dth * dph));
    eta_t = (16/pi) * fD^2 * abs(I_num)^2 / I_den;

    D_dBi(i_fD) = 10 * log10(eta_t * D_max);
    G_dBi(i_fD) = 10 * log10(eta_s * eta_t * D_max);
end

%% Plot vs D
D_range = f_focal ./ fD_range;

fig = figure('Color', 'w', 'Position', [100 100 800 550]);
set(fig, 'InvertHardcopy', 'off');
plot(D_range, D_max_dBi, 'b-', 'LineWidth', 1.5); hold on;
plot(D_range, D_dBi, 'r--', 'LineWidth', 1.5);
plot(D_range, G_dBi, 'k-', 'LineWidth', 2);
xlabel('D [m]');
ylabel('[dBi]');
title('Q5: Directivity and Gain vs D (f = 5 m, 500 GHz)');
legend('D_{max} (uniform aperture)', 'Directivity D', 'Gain G', 'Location', 'best');
grid on;
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q5_directivity_gain.png'));
