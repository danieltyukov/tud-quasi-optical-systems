%% Q4: Efficiencies vs Diameter (3 points)

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

%% Pre-compute feed pattern on 2D grid
Nth_feed = 901;
Nph_feed = 361;
theta_feed = linspace(0, pi/2, Nth_feed);
phi_feed = linspace(0, 2*pi, Nph_feed);
dth = theta_feed(2) - theta_feed(1);
dph = phi_feed(2) - phi_feed(1);
[TH_feed, PH_feed] = meshgrid(theta_feed, phi_feed);

[E_th_2d, E_ph_2d] = feed_farfield(k0, a_feed, TH_feed, PH_feed);

U_feed = abs(E_th_2d).^2 + abs(E_ph_2d).^2;
E_co = E_th_2d .* sin(PH_feed) + E_ph_2d .* cos(PH_feed); % Ludwig-3 co-pol
G_feed = U_feed;
P_total = sum(sum(U_feed .* sin(TH_feed) * dth * dph));

%% Sweep f/D
N_fD = 200;
fD_range = linspace(0.6, 6, N_fD);
eta_s = zeros(1, N_fD);
eta_t = zeros(1, N_fD);
eta_ap = zeros(1, N_fD);

for i_fD = 1:N_fD
    fD = fD_range(i_fD);
    D_refl = f_focal / fD;
    theta_0 = 2 * atan(D_refl / (4 * f_focal));
    idx_th0 = find(theta_feed <= theta_0, 1, 'last');

    P_intercepted = sum(sum(U_feed(:, 1:idx_th0) .* sin(TH_feed(:, 1:idx_th0)) * dth * dph));
    eta_s(i_fD) = P_intercepted / P_total;

    integrand_num = abs(E_co(:, 1:idx_th0)) .* tan(TH_feed(:, 1:idx_th0)/2);
    I_num = sum(sum(integrand_num * dth * dph));
    integrand_den = G_feed(:, 1:idx_th0) .* sin(TH_feed(:, 1:idx_th0));
    I_den = sum(sum(integrand_den * dth * dph));
    eta_t(i_fD) = (16/pi) * fD^2 * abs(I_num)^2 / I_den;

    eta_ap(i_fD) = eta_s(i_fD) * eta_t(i_fD);
end

%% Plot vs D
D_range = f_focal ./ fD_range;
[eta_ap_max, idx_max] = max(eta_ap);
fD_opt = fD_range(idx_max);
D_opt = f_focal / fD_opt;

fig = figure('Color', 'w', 'Position', [100 100 800 550]);
set(fig, 'InvertHardcopy', 'off');
plot(D_range, eta_s, 'b-', 'LineWidth', 1.5); hold on;
plot(D_range, eta_t, 'r--', 'LineWidth', 1.5);
plot(D_range, eta_ap, 'k-', 'LineWidth', 2);
plot(D_opt, eta_ap_max, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(D_opt + 0.1, eta_ap_max - 0.05, sprintf('max \\eta_{ap} = %.2f\nf/D = %.2f, D = %.2f m', ...
     eta_ap_max, fD_opt, D_opt), 'FontSize', 10);
xlabel('D [m]');
ylabel('Efficiency');
title('Q4: Efficiencies vs D (f = 5 m, 500 GHz)');
legend('Spillover \eta_s', 'Taper \eta_t', 'Aperture \eta_{ap}', 'Location', 'best');
grid on;
ylim([0 1]);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q4_efficiencies.png'));

fprintf('Max aperture efficiency: %.4f at f/D = %.2f\n', eta_ap_max, fD_opt);
