%% Q1.3-Q1.6: Full Rx analysis — Voc, pattern, Tx comparison
clear; close all; clc;
set_plot_defaults();

%% Parameters
c = 3e8;
freq = 220e9;
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
zeta0 = 120 * pi;
E0_PW = 1;

D_f = 3.6 * lambda0;
a_feed = D_f / 2;
D_r = 130 * lambda0;
a_refl = D_r / 2;
f_num = 1.8;
F = f_num * D_r;
theta0 = 2 * atan(D_r / (4 * F));
A_r = pi * a_refl^2;

fprintf('=== System Parameters ===\n');
fprintf('f = %.0f GHz, lambda0 = %.4f mm\n', freq/1e9, lambda0*1e3);
fprintf('D_f = %.2f lambda, D_r = %.0f lambda\n', D_f/lambda0, D_r/lambda0);
fprintf('F = %.4f mm, f# = %.1f\n', F*1e3, f_num);
fprintf('theta0 = %.2f deg\n', rad2deg(theta0));

%% Pre-compute feed far field on 2D grid over FO sphere
Nth = 401;
Nph = 361;
theta_arr = linspace(0, pi/2, Nth);
phi_arr = linspace(0, 2*pi, Nph);
dth = theta_arr(2) - theta_arr(1);
dph = phi_arr(2) - phi_arr(1);
[TH, PH] = meshgrid(theta_arr, phi_arr);

[Eth_feed, Eph_feed] = feed_farfield_x(k0, a_feed, TH, PH);

% Total radiated power of feed (integrating over upper hemisphere)
U_feed = abs(Eth_feed).^2 + abs(Eph_feed).^2;
P_rad = (1/(2*zeta0)) * sum(U_feed .* sin(TH) * dth * dph, 'all');

fprintf('P_rad (feed) = %.6e W\n', P_rad);

%% Q1.3: Broadside reception — Voc and aperture efficiency

% GO field on FO sphere (broadside, x-pol PW)
idx_th0 = find(theta_arr <= theta0, 1, 'last');

TH_fo = TH(:, 1:idx_th0);
PH_fo = PH(:, 1:idx_th0);
Eth_f = Eth_feed(:, 1:idx_th0);
Eph_f = Eph_feed(:, 1:idx_th0);

S_par = 2 ./ (1 + cos(TH_fo));
Eth_GO = -S_par .* cos(PH_fo) * E0_PW;
Eph_GO = S_par .* sin(PH_fo) * E0_PW;

% Voc = (2F/zeta0) * integral of (e_a . e_GO) sin(theta) dtheta dphi
integrand_Voc = (Eth_f .* Eth_GO + Eph_f .* Eph_GO) .* sin(TH_fo);
Voc = (2 * F / zeta0) * sum(integrand_Voc * dth * dph, 'all');

P_load = abs(Voc)^2 / (16 * P_rad);
P_inc = (abs(E0_PW)^2 / (2 * zeta0)) * A_r;
eta_ap_rx = P_load / P_inc;

fprintf('\n=== Q1.3: Broadside Reception ===\n');
fprintf('|Voc| = %.6e V\n', abs(Voc));
fprintf('P_load = %.6e W\n', P_load);
fprintf('P_inc = %.6e W\n', P_inc);
fprintf('eta_ap (Rx) = %.4f (%.2f%%)\n', eta_ap_rx, eta_ap_rx*100);

% Decompose into spillover and taper
P_feed_in_cone = (1/(2*zeta0)) * sum(U_feed(:, 1:idx_th0) .* sin(TH_fo) * dth * dph, 'all');
eta_so_rx = P_feed_in_cone / P_rad;

% Taper efficiency from the Cauchy-Schwarz perspective
num_sq = abs(sum(integrand_Voc * dth * dph, 'all'))^2;
denom1 = sum(U_feed(:, 1:idx_th0) .* sin(TH_fo) * dth * dph, 'all');
GO_power = abs(Eth_GO).^2 + abs(Eph_GO).^2;
denom2 = sum(GO_power .* sin(TH_fo) * dth * dph, 'all');
eta_t_rx = num_sq / (denom1 * denom2);

fprintf('eta_so (Rx) = %.4f (%.2f%%)\n', eta_so_rx, eta_so_rx*100);
fprintf('eta_t (Rx) = %.4f (%.2f%%)\n', eta_t_rx, eta_t_rx*100);
fprintf('eta_so * eta_t = %.4f\n', eta_so_rx * eta_t_rx);

%% Q1.4 + Q1.5: Rx radiation pattern vs theta_i

theta_i_max = 3; % degrees
Nth_i = 301;
theta_i_arr = linspace(0, deg2rad(theta_i_max), Nth_i);

% Precompute rho_refl = 2F*tan(theta/2) for the off-broadside phase
rho_refl = 2 * F * tan(TH_fo / 2);

% Broadside dot product (without off-broadside phase) for reference
dot_bs = Eth_f .* Eth_GO + Eph_f .* Eph_GO;

% Rx pattern in E-plane (phi_i = 0) and H-plane (phi_i = 90)
P_rx_E = zeros(1, Nth_i);
P_rx_H = zeros(1, Nth_i);

for ii = 1:Nth_i
    thi = theta_i_arr(ii);
    dk_rho = k0 * sin(thi);

    % E-plane: phi_i = 0
    phase_E = -dk_rho * rho_refl .* cos(PH_fo - 0);
    integrand_E = dot_bs .* exp(1j * phase_E) .* sin(TH_fo);
    Voc_E = sum(integrand_E * dth * dph, 'all');
    P_rx_E(ii) = abs(Voc_E)^2;

    % H-plane: phi_i = pi/2
    phase_H = -dk_rho * rho_refl .* cos(PH_fo - pi/2);
    integrand_H = dot_bs .* exp(1j * phase_H) .* sin(TH_fo);
    Voc_H = sum(integrand_H * dth * dph, 'all');
    P_rx_H(ii) = abs(Voc_H)^2;
end

P_rx_E_dB = 10 * log10(P_rx_E / P_rx_E(1));
P_rx_H_dB = 10 * log10(P_rx_H / P_rx_H(1));

%% Tx far field for comparison
fprintf('\nComputing Tx far field for comparison...\n');

N_ap = 261;
x_ap = linspace(-a_refl, a_refl, N_ap);
y_ap = linspace(-a_refl, a_refl, N_ap);
dx_ap = x_ap(2) - x_ap(1);
dy_ap = y_ap(2) - y_ap(1);
[X_ap, Y_ap] = meshgrid(x_ap, y_ap);
rho_ap = sqrt(X_ap.^2 + Y_ap.^2);
mask = rho_ap <= a_refl;

theta_p = 2 * atan(rho_ap / (2 * F));
phi_p = atan2(Y_ap, X_ap);
phi_p(rho_ap < 1e-15) = 0;

% Aperture field from x-polarized feed
[Eth_f_ap, Eph_f_ap] = feed_farfield_x(k0, a_feed, theta_p, phi_p);
amp_taper = (1 + cos(theta_p)) ./ (2 * F);

% Project reflected field onto aperture plane (x,y components)
E_a_x = amp_taper .* (Eth_f_ap .* cos(theta_p) .* cos(phi_p) - Eph_f_ap .* sin(phi_p));
E_a_y = amp_taper .* (Eth_f_ap .* cos(theta_p) .* sin(phi_p) + Eph_f_ap .* cos(phi_p));
E_a_x = E_a_x .* mask;
E_a_y = E_a_y .* mask;

% Equivalent currents
M_x = E_a_y;
M_y = -E_a_x;
J_x = -E_a_x / zeta0;
J_y = -E_a_y / zeta0;

% Tx far field at theta_i angles
E_tx_E = zeros(1, Nth_i);
E_tx_H = zeros(1, Nth_i);

for ii = 1:Nth_i
    th = theta_i_arr(ii);

    % E-plane (phi_obs = 0 for x-pol)
    kx_ff = k0 * sin(th);
    ky_ff = 0;
    phase_ap = exp(1j * (kx_ff * X_ap + ky_ff * Y_ap));

    Nx = sum(J_x .* phase_ap, 'all') * dx_ap * dy_ap;
    Ny = sum(J_y .* phase_ap, 'all') * dx_ap * dy_ap;
    Lx = sum(M_x .* phase_ap, 'all') * dx_ap * dy_ap;
    Ly = sum(M_y .* phase_ap, 'all') * dx_ap * dy_ap;

    cos_th = cos(th);
    E_th = Ly + zeta0 * Nx * cos_th;  % phi=0: sin_ph=0, cos_ph=1
    E_ph = Lx * cos_th - zeta0 * Ny;
    E_tx_E(ii) = sqrt(abs(E_th)^2 + abs(E_ph)^2);

    % H-plane (phi_obs = pi/2 for x-pol)
    kx_ff = 0;
    ky_ff = k0 * sin(th);
    phase_ap = exp(1j * (kx_ff * X_ap + ky_ff * Y_ap));

    Nx = sum(J_x .* phase_ap, 'all') * dx_ap * dy_ap;
    Ny = sum(J_y .* phase_ap, 'all') * dx_ap * dy_ap;
    Lx = sum(M_x .* phase_ap, 'all') * dx_ap * dy_ap;
    Ly = sum(M_y .* phase_ap, 'all') * dx_ap * dy_ap;

    cos_th = cos(th);
    E_th = -Lx + zeta0 * Ny * cos_th; % phi=pi/2: sin_ph=1, cos_ph=0
    E_ph = Ly * cos_th + zeta0 * Nx;
    E_tx_H(ii) = sqrt(abs(E_th)^2 + abs(E_ph)^2);
end

E_tx_E_dB = 20 * log10(E_tx_E / max(E_tx_E));
E_tx_H_dB = 20 * log10(E_tx_H / max(E_tx_H));

%% Plot Rx vs Tx pattern
fig = figure('Color', 'w', 'Position', [100 100 1000 450]);
set(fig, 'InvertHardcopy', 'off');

subplot(1,2,1);
plot(rad2deg(theta_i_arr), P_rx_E_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(rad2deg(theta_i_arr), E_tx_E_dB, 'r--', 'LineWidth', 1.5);
xlabel('\theta_i [deg]');
ylabel('Normalized pattern [dB]');
title('\phi = 0° (E-plane)');
legend('Rx analysis', 'Tx analysis', 'Location', 'northeast');
grid on; ylim([-40 0]); xlim([0 3]);
set(gca, 'Color', 'w');

subplot(1,2,2);
plot(rad2deg(theta_i_arr), P_rx_H_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(rad2deg(theta_i_arr), E_tx_H_dB, 'r--', 'LineWidth', 1.5);
xlabel('\theta_i [deg]');
ylabel('Normalized pattern [dB]');
title('\phi = 90° (H-plane)');
legend('Rx analysis', 'Tx analysis', 'Location', 'northeast');
grid on; ylim([-40 0]); xlim([0 3]);
set(gca, 'Color', 'w');

sgtitle('Q1.5: Rx vs Tx Radiation Pattern (220 GHz, f_# = 1.8)');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_rx_vs_tx.png'));

%% Q1.6: Tx/Rx comparison of efficiency, directivity, gain

% Tx efficiencies
eta_so_tx = P_feed_in_cone / P_rad;  % same as Rx

% Taper: using Ludwig-3 co-pol for x-polarized feed
E_co_tx = Eth_feed(:, 1:idx_th0) .* cos(PH_fo) - Eph_feed(:, 1:idx_th0) .* sin(PH_fo);
I_num_tx = sum(abs(E_co_tx) .* tan(TH_fo/2) * dth * dph, 'all');
I_den_tx = sum(U_feed(:, 1:idx_th0) .* sin(TH_fo) * dth * dph, 'all');
eta_t_tx = (16/pi) * f_num^2 * abs(I_num_tx)^2 / I_den_tx;

eta_ap_tx = eta_so_tx * eta_t_tx;

% For Rx, define taper consistently: eta_t_rx = eta_ap_rx / eta_so_rx
eta_t_rx_consistent = eta_ap_rx / eta_so_rx;

D_max = (pi * D_r / lambda0)^2;
D_max_dBi = 10 * log10(D_max);

% Directivity = D_max * eta_t, Gain = D_max * eta_ap
D_tx = eta_t_tx * D_max;
D_rx = eta_t_rx_consistent * D_max;
G_tx = eta_ap_tx * D_max;
G_rx = eta_ap_rx * D_max;

fprintf('\n=== Q1.6: Tx vs Rx Comparison ===\n');
fprintf('                    Tx          Rx\n');
fprintf('eta_spillover:   %.4f       %.4f\n', eta_so_tx, eta_so_rx);
fprintf('eta_taper:       %.4f       %.4f\n', eta_t_tx, eta_t_rx_consistent);
fprintf('eta_aperture:    %.4f       %.4f\n', eta_ap_tx, eta_ap_rx);
fprintf('D_max:           %.2f dBi\n', D_max_dBi);
fprintf('Directivity:     %.2f dBi   %.2f dBi\n', 10*log10(D_tx), 10*log10(D_rx));
fprintf('Gain:            %.2f dBi   %.2f dBi\n', 10*log10(G_tx), 10*log10(G_rx));

% Save workspace for Q3 displaced feed analysis
save('rx_workspace.mat', 'k0', 'lambda0', 'zeta0', 'a_feed', 'D_r', 'a_refl', ...
     'F', 'f_num', 'theta0', 'A_r', 'E0_PW', ...
     'theta_arr', 'phi_arr', 'dth', 'dph', 'TH', 'PH', ...
     'Eth_feed', 'Eph_feed', 'P_rad', 'U_feed', ...
     'idx_th0', 'TH_fo', 'PH_fo', 'Eth_f', 'Eph_f', ...
     'S_par', 'Eth_GO', 'Eph_GO', 'dot_bs', 'rho_refl', ...
     'P_rx_E', 'theta_i_arr', 'Nth_i', 'eta_ap_rx');
