%% Assignment 4 - Part 1: Feed at the center of the focal plane
%  FO for FPA - Parabolic Reflector Analysis in Reception
%
%  Circular feed with x-polarized uniform amplitude, centered at focal plane.

clear; close all; clc;
set_plot_defaults();

%% Constants and geometry
c = 3e8;
freq = 220e9;                       % 220 GHz
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
zeta = 120 * pi;                    % free-space impedance

D_f = 3.6 * lambda0;               % feed diameter
a_feed = D_f / 2;                   % feed radius

D_r = 130 * lambda0;               % reflector diameter
a_refl = D_r / 2;
f_num = 1.8;                       % f-number = F/D
F = f_num * D_r;                   % focal length

theta_0 = 2 * atan(D_r / (4*F));   % reflector rim angle
% Also: theta_0 = 2*atan(1/(4*f_num))

E0_PW = 1;                         % incident PW amplitude [V/m]

fprintf('=== Assignment 4 - Part 1: FO for FPA ===\n');
fprintf('Frequency: %.0f GHz, lambda0 = %.4f mm\n', freq/1e9, lambda0*1e3);
fprintf('Feed: D_f = %.1f lambda0 = %.4f mm\n', D_f/lambda0, D_f*1e3);
fprintf('Reflector: D_r = %.0f lambda0 = %.4f m\n', D_r/lambda0, D_r);
fprintf('f/# = %.1f, F = %.4f m\n', f_num, F);
fprintf('theta_0 = %.4f deg\n', theta_0*180/pi);
fprintf('HPBW (Airy): Delta_theta = lambda/D = %.5f deg\n', (lambda0/D_r)*180/pi);
fprintf('First null: 1.22*lambda/D = %.5f deg\n', 1.22*(lambda0/D_r)*180/pi);

%% Angular grid on FO sphere
Nth = 501;
Nph = 361;
theta_vec = linspace(0, pi/2, Nth);
phi_vec = linspace(0, 2*pi, Nph);
dth = theta_vec(2) - theta_vec(1);
dph = phi_vec(2) - phi_vec(1);
[TH, PH] = meshgrid(theta_vec, phi_vec);  % size: Nph x Nth

%% ========== Q1.1: GO electric field on FO sphere (broadside) ==========
fprintf('\n--- Q1.1: GO field on FO sphere ---\n');

[e_GO_th, e_GO_ph] = GO_field_parabolic(E0_PW, TH, PH, theta_0);

% The "voltage" pattern on the FO sphere: V_r^GO = e_r^GO * R * exp(jkR)
% For field matching, we use R=F as the FO sphere radius.
% The R*exp(jkR) factor will cancel in the Voc computation, so we can
% work with just the angular field pattern e_r^GO directly as V_r^GO.
% (We set R=1 for normalization consistency.)
VGO_th = e_GO_th;
VGO_ph = e_GO_ph;

fprintf('GO field computed on %d x %d grid\n', Nph, Nth);
fprintf('Max |e_GO_theta| = %.4f V/m\n', max(abs(e_GO_th(:))));
fprintf('Max |e_GO_phi| = %.4f V/m\n', max(abs(e_GO_ph(:))));

%% ========== Q1.2: Feed far field on FO sphere & comparison ==========
fprintf('\n--- Q1.2: Feed far field vs GO field ---\n');

% Feed far field (x-polarized circular aperture)
[E_feed_th, E_feed_ph] = feed_farfield_x(k0, a_feed, TH, PH);

% The feed "voltage" pattern: V_a^tx = E_a^tx * r * exp(jkr)
% Again R_FF=1 in feed_farfield_x, so the output is already V_a^tx.
Va_th = E_feed_th;
Va_ph = E_feed_ph;

% --- Plot phi=0 plane ---
idx_phi0 = 1;  % phi=0
GO_total_phi0 = sqrt(abs(VGO_th(idx_phi0,:)).^2 + abs(VGO_ph(idx_phi0,:)).^2);
Feed_total_phi0 = sqrt(abs(Va_th(idx_phi0,:)).^2 + abs(Va_ph(idx_phi0,:)).^2);

% Normalize to their own maxima for comparison
GO_max = max(GO_total_phi0);
Feed_max = max(Feed_total_phi0);
GO_phi0_dB = 20*log10(GO_total_phi0 / GO_max + 1e-30);
Feed_phi0_dB = 20*log10(Feed_total_phi0 / Feed_max + 1e-30);

fig1 = figure('Color', 'w', 'Position', [100 100 900 500]);
subplot(1,2,1);
plot(theta_vec*180/pi, GO_phi0_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_vec*180/pi, Feed_phi0_dB, 'r--', 'LineWidth', 1.5);
xline(theta_0*180/pi, 'k:', sprintf('\\theta_0 = %.1f°', theta_0*180/pi), 'LineWidth', 1);
xlabel('\theta [deg]');
ylabel('Normalized field [dB]');
title('\phi = 0° plane (H-plane)');
legend('GO field', 'Feed far field', 'Location', 'southwest');
grid on; ylim([-40 0]); xlim([0 45]);

% --- Plot phi=90 plane ---
idx_phi90 = find(abs(phi_vec - pi/2) < dph/2, 1);
GO_total_phi90 = sqrt(abs(VGO_th(idx_phi90,:)).^2 + abs(VGO_ph(idx_phi90,:)).^2);
Feed_total_phi90 = sqrt(abs(Va_th(idx_phi90,:)).^2 + abs(Va_ph(idx_phi90,:)).^2);

GO_max90 = max(GO_total_phi90);
Feed_max90 = max(Feed_total_phi90);
GO_phi90_dB = 20*log10(GO_total_phi90 / GO_max90 + 1e-30);
Feed_phi90_dB = 20*log10(Feed_total_phi90 / Feed_max90 + 1e-30);

subplot(1,2,2);
plot(theta_vec*180/pi, GO_phi90_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_vec*180/pi, Feed_phi90_dB, 'r--', 'LineWidth', 1.5);
xline(theta_0*180/pi, 'k:', sprintf('\\theta_0 = %.1f°', theta_0*180/pi), 'LineWidth', 1);
xlabel('\theta [deg]');
ylabel('Normalized field [dB]');
title('\phi = 90° plane (E-plane)');
legend('GO field', 'Feed far field', 'Location', 'southwest');
grid on; ylim([-40 0]); xlim([0 45]);

sgtitle('Q1.2: GO Field vs Feed Far Field on FO Sphere');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_2_GO_vs_feed.png'));

%% ========== Q1.3: Power received & aperture efficiency (broadside) ==========
fprintf('\n--- Q1.3: Power received and aperture efficiency ---\n');

% Open-circuit voltage via field matching (reaction integral)
Voc = compute_Voc(Va_th, Va_ph, VGO_th, VGO_ph, TH, dth, dph, zeta);

% Radiated power of the feed antenna (for normalization)
% P_rad = (1/zeta) * integral{ |V_a^tx|^2 * sin(theta) dtheta dphi }
% Using the pattern over the full forward hemisphere
U_feed = abs(Va_th).^2 + abs(Va_ph).^2;
P_rad_feed = (1/zeta) * sum(U_feed(:) .* sin(TH(:))) * dth * dph;

% Also compute over full sphere (add back hemisphere contribution)
% For a circular aperture the back radiation is negligible, but for
% completeness we integrate over [0, pi/2] and double (symmetry assumption).
% Actually, the feed radiates into forward hemisphere only (aperture antenna).
% P_rad = integral over [0, pi] but the aperture pattern is zero for theta > pi/2.
% So P_rad_feed computed above is the total radiated power.

% Power delivered to load (matched: chi_match = 1)
% P_load = |Voc|^2 / (8*R_a) where R_a = 2*P_rad / |I0|^2, with I0=1
% So P_load = |Voc|^2 / (16*P_rad)
P_load = abs(Voc)^2 / (16 * P_rad_feed);

% Incident power on the reflector
A_r = pi * a_refl^2;
P_inc = (abs(E0_PW)^2 / (2*zeta)) * A_r;

% Aperture efficiency
eta_ap = P_load / P_inc;

fprintf('|Voc| = %.6e\n', abs(Voc));
fprintf('P_rad (feed) = %.6e W\n', P_rad_feed);
fprintf('P_load = %.6e W\n', P_load);
fprintf('P_inc = %.6e W\n', P_inc);
fprintf('Aperture efficiency eta_ap = %.4f (%.2f%%)\n', eta_ap, eta_ap*100);

% Decompose efficiency: eta_ap = eta_so * eta_t * chi_match
% Spillover efficiency: fraction of feed Tx power intercepted by reflector
idx_th0 = find(theta_vec <= theta_0, 1, 'last');
P_intercepted = (1/zeta) * sum(sum(U_feed(:,1:idx_th0) .* sin(TH(:,1:idx_th0)))) * dth * dph;
eta_so = P_intercepted / P_rad_feed;

% Taper efficiency (from Voc formula)
eta_t = eta_ap / eta_so;  % since chi_match = 1

fprintf('Spillover efficiency eta_so = %.4f\n', eta_so);
fprintf('Taper efficiency eta_t = %.4f\n', eta_t);

%% ========== Q1.4 & Q1.5: Scan over incident angles ==========
fprintf('\n--- Q1.4/Q1.5: Scanning over incident angles ---\n');

% Incident angle range
N_thi = 61;
N_phi_i = 73;
theta_i_vec = linspace(0, 3*pi/180, N_thi);    % 0 to 3 degrees
phi_i_vec = linspace(0, 2*pi, N_phi_i);         % 0 to 360 degrees
dth_i = theta_i_vec(2) - theta_i_vec(1);

P_rx = zeros(N_thi, N_phi_i);

fprintf('Computing Rx pattern for %d x %d incident angles...\n', N_thi, N_phi_i);

for i_thi = 1:N_thi
    for i_phi = 1:N_phi_i
        thi = theta_i_vec(i_thi);
        phi_i = phi_i_vec(i_phi);

        % GO field for this incidence angle
        [eGO_th_off, eGO_ph_off] = GO_field_off_broadside(E0_PW, k0, F, ...
                                        TH, PH, theta_0, thi, phi_i);

        % Open-circuit voltage
        Voc_ij = compute_Voc(Va_th, Va_ph, eGO_th_off, eGO_ph_off, TH, dth, dph, zeta);

        % Received power (matched)
        P_rx(i_thi, i_phi) = abs(Voc_ij)^2 / (16 * P_rad_feed);
    end

    if mod(i_thi, 10) == 0
        fprintf('  theta_i = %.2f deg done\n', thi*180/pi);
    end
end

% Normalize and convert to dB
P_rx_max = max(P_rx(:));
P_rx_norm_dB = 10 * log10(P_rx / P_rx_max + 1e-30);

% Plot: Rx pattern in phi_i=0 plane (1D cut)
idx_phi0_i = 1;
idx_phi180_i = find(abs(phi_i_vec - pi) < (phi_i_vec(2)-phi_i_vec(1))/2, 1);

fig2 = figure('Color', 'w', 'Position', [100 100 900 500]);
plot(theta_i_vec*180/pi, P_rx_norm_dB(:, idx_phi0_i), 'b-', 'LineWidth', 1.5); hold on;

% Also plot phi_i = 90 deg cut
idx_phi90_i = find(abs(phi_i_vec - pi/2) < (phi_i_vec(2)-phi_i_vec(1))/2, 1);
if ~isempty(idx_phi90_i)
    plot(theta_i_vec*180/pi, P_rx_norm_dB(:, idx_phi90_i), 'r--', 'LineWidth', 1.5);
end

% Reference: Airy pattern for uniform circular aperture
airy_ref = ones(size(theta_i_vec));
for it = 1:N_thi
    arg = k0 * a_refl * sin(theta_i_vec(it));
    if abs(arg) > 1e-10
        airy_ref(it) = abs(2 * besselj(1, arg) / arg);
    end
end
airy_dB = 20*log10(airy_ref / max(airy_ref));
plot(theta_i_vec*180/pi, airy_dB, 'k:', 'LineWidth', 1.5);

xlabel('\theta_i [deg]');
ylabel('Normalized received power [dB]');
title('Q1.5: Rx Pattern of Reflector Antenna (centered feed)');
legend('\phi_i = 0° (H-plane)', '\phi_i = 90° (E-plane)', 'Airy pattern (ref)', ...
       'Location', 'northeast');
grid on; ylim([-40 0]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_5_rx_pattern.png'));

%% ========== Q1.6: Compare Tx vs Rx aperture efficiencies ==========
fprintf('\n--- Q1.6: Tx vs Rx comparison ---\n');

% Maximum directivity (uniform aperture)
D_max = (pi * D_r / lambda0)^2;
D_max_dBi = 10*log10(D_max);

% Rx analysis results
D_rx_dBi = 10*log10(eta_t * D_max);
G_rx_dBi = 10*log10(eta_ap * D_max);

fprintf('--- Rx Analysis ---\n');
fprintf('D_max = %.2f dBi\n', D_max_dBi);
fprintf('eta_ap (Rx) = %.4f\n', eta_ap);
fprintf('eta_so (Rx) = %.4f\n', eta_so);
fprintf('eta_t  (Rx) = %.4f\n', eta_t);
fprintf('Directivity (Rx) = %.2f dBi\n', D_rx_dBi);
fprintf('Gain (Rx) = %.2f dBi\n', G_rx_dBi);

% Tx analysis (from Assignment 2 approach - compute on the fly)
% Taper efficiency in Tx formulation:
% eta_t = (1/A) * |integral E_a^co dA|^2 / integral |E_a|^2 dA
% Using the angular formulation with cos^2(theta'/2) feed pattern:
% eta_t_tx = 32*(F/D)^2 * |integral sqrt(G_co)*tan(theta/2) dtheta|^2 /
%                          integral G(theta)*sin(theta) dtheta

% Compute feed radiated power pattern for Tx analysis
% For the x-polarized feed, the co-pol component depends on cut:
% Ludwig-3 co-pol for x-polarized: E_co = E_th*cos(phi) - E_ph*sin(phi)
E_co_feed = Va_th .* cos(PH) - Va_ph .* sin(PH);
G_feed = U_feed;

% Spillover (Tx) - same as Rx by reciprocity
eta_so_tx = eta_so;

% Taper (Tx)
integrand_num_tx = abs(E_co_feed(:,1:idx_th0)) .* tan(TH(:,1:idx_th0)/2);
I_num_tx = sum(integrand_num_tx(:)) * dth * dph;
integrand_den_tx = G_feed(:,1:idx_th0) .* sin(TH(:,1:idx_th0));
I_den_tx = sum(integrand_den_tx(:)) * dth * dph;
eta_t_tx = 32 * f_num^2 * abs(I_num_tx)^2 / I_den_tx;

eta_ap_tx = eta_so_tx * eta_t_tx;
D_tx_dBi = 10*log10(eta_t_tx * D_max);
G_tx_dBi = 10*log10(eta_ap_tx * D_max);

fprintf('\n--- Tx Analysis ---\n');
fprintf('eta_ap (Tx) = %.4f\n', eta_ap_tx);
fprintf('eta_so (Tx) = %.4f\n', eta_so_tx);
fprintf('eta_t  (Tx) = %.4f\n', eta_t_tx);
fprintf('Directivity (Tx) = %.2f dBi\n', D_tx_dBi);
fprintf('Gain (Tx) = %.2f dBi\n', G_tx_dBi);

fprintf('\n--- Comparison ---\n');
fprintf('Delta eta_ap = %.6f\n', abs(eta_ap - eta_ap_tx));
fprintf('Delta Gain = %.4f dB\n', abs(G_rx_dBi - G_tx_dBi));
fprintf('(Should be identical by reciprocity)\n');

%% Save workspace for Part 2
save(fullfile('..', 'matlab', 'Part1_workspace.mat'), ...
     'lambda0', 'k0', 'zeta', 'D_f', 'a_feed', 'D_r', 'a_refl', ...
     'f_num', 'F', 'theta_0', 'E0_PW', 'freq', ...
     'theta_vec', 'phi_vec', 'dth', 'dph', 'TH', 'PH', ...
     'Va_th', 'Va_ph', 'P_rad_feed', 'eta_ap', 'eta_so', 'eta_t', ...
     'P_rx', 'theta_i_vec', 'phi_i_vec', 'P_rx_max');

fprintf('\nPart 1 workspace saved.\n');
