%% Q2.1 + Q2.2: Displaced feed analysis
clear; close all; clc;
set_plot_defaults();

%% Load precomputed data from Q2
load('rx_workspace.mat');

D_f = 3.6 * lambda0;

%% Q2.1: Displaced feed, dx = Df = 3.6*lambda0

dx_disp = D_f;  % = 3.6*lambda0 = 2*lambda0*f_num (since f_num=1.8)

% Displaced feed adds phase: exp(j*k0*dx*sin(theta)*cos(phi))
disp_phase = exp(1j * k0 * dx_disp * sin(TH_fo) .* cos(PH_fo));

dot_disp = (Eth_f .* disp_phase .* Eth_GO + Eph_f .* disp_phase .* Eph_GO);

% Rx pattern for displaced feed
P_rx_disp_E = zeros(1, Nth_i);
P_rx_cen_E = zeros(1, Nth_i);

theta_i_max_scan = 3;
theta_i_scan = linspace(-deg2rad(theta_i_max_scan), deg2rad(theta_i_max_scan), 2*Nth_i-1);

P_rx_disp_scan = zeros(size(theta_i_scan));
P_rx_cen_scan = zeros(size(theta_i_scan));

for ii = 1:length(theta_i_scan)
    thi = theta_i_scan(ii);
    dk_rho = k0 * sin(abs(thi));
    phi_i = 0;
    if thi < 0
        phi_i = pi;
        dk_rho = k0 * sin(abs(thi));
    end

    phase_OB = -dk_rho * rho_refl .* cos(PH_fo - phi_i);

    % Central feed
    integrand_cen = dot_bs .* exp(1j * phase_OB) .* sin(TH_fo);
    P_rx_cen_scan(ii) = abs(sum(integrand_cen * dth * dph, 'all'))^2;

    % Displaced feed
    integrand_disp = dot_disp .* exp(1j * phase_OB) .* sin(TH_fo);
    P_rx_disp_scan(ii) = abs(sum(integrand_disp * dth * dph, 'all'))^2;
end

P_max_cen = max(P_rx_cen_scan);
P_rx_cen_dB = 10 * log10(P_rx_cen_scan / P_max_cen);
P_rx_disp_dB = 10 * log10(P_rx_disp_scan / P_max_cen);

% Find pointing direction of displaced feed
[P_disp_max, idx_disp_max] = max(P_rx_disp_scan);
theta_point = rad2deg(theta_i_scan(idx_disp_max));

% Aperture efficiency of displaced feed
eta_ap_disp = (P_disp_max / P_max_cen) * eta_ap_rx;

fprintf('=== Q2.1: Displaced Feed (dx = %.1f lambda0 = %.1f mm) ===\n', dx_disp/lambda0, dx_disp*1e3);
fprintf('Pointing direction: theta_s = %.3f deg\n', theta_point);
fprintf('Theoretical: theta_s = arcsin(dx/F) = %.3f deg\n', rad2deg(asin(dx_disp/F)));
fprintf('Aperture efficiency (displaced): %.4f (%.2f%%)\n', eta_ap_disp, eta_ap_disp*100);
fprintf('Aperture efficiency (central):   %.4f (%.2f%%)\n', eta_ap_rx, eta_ap_rx*100);
fprintf('Scan loss: %.2f dB\n', 10*log10(P_disp_max/P_max_cen));

%% Plot Q2.1
fig1 = figure('Color', 'w', 'Position', [100 100 800 500]);
set(fig1, 'InvertHardcopy', 'off');

plot(rad2deg(theta_i_scan), P_rx_cen_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(rad2deg(theta_i_scan), P_rx_disp_dB, 'r--', 'LineWidth', 1.5);
xlabel('\theta_i [deg]');
ylabel('Normalized P_{rx} [dB]');
title(sprintf('Q2.1: Central vs Displaced Feed (d_x = %.1f\\lambda_0)', dx_disp/lambda0));
legend('Central feed', sprintf('Displaced (d_x = %.1f\\lambda_0)', dx_disp/lambda0), ...
       'Location', 'best');
grid on; ylim([-40 0]);
xlim([-theta_i_max_scan theta_i_max_scan]);
set(gca, 'Color', 'w');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q3_displaced_feed.png'));

%% Q2.2: Three sampling schemes

dx_cases = [2*lambda0*f_num, lambda0*f_num, 0.5*lambda0*f_num];
case_labels = {'2\lambda_0 f_#', '\lambda_0 f_#', '0.5\lambda_0 f_#'};
theoretical_cross = [-3, -3, -3]; % Airy beam crossing for each spacing
% Actually: 2*lambda0*f# -> field sampling -> -3 dB crossing
%           lambda0*f# -> HPBW -> depends
%           0.5*lambda0*f# -> power sampling -> -1.5 dB crossing
% But the crossing depends on where beams overlap

% Extended angular scan for sampling analysis
theta_i_ext = linspace(-deg2rad(4), deg2rad(4), 801);

fig2 = figure('Color', 'w', 'Position', [100 100 1200 400]);
set(fig2, 'InvertHardcopy', 'off');
colors = {'r', [0 0.6 0], [0.8 0.4 0]};

for ic = 1:3
    dx = dx_cases(ic);
    disp_ph = exp(1j * k0 * dx * sin(TH_fo) .* cos(PH_fo));
    dot_d = (Eth_f .* disp_ph .* Eth_GO + Eph_f .* disp_ph .* Eph_GO);

    P_cen = zeros(size(theta_i_ext));
    P_disp = zeros(size(theta_i_ext));

    for ii = 1:length(theta_i_ext)
        thi = theta_i_ext(ii);
        phi_i = 0;
        if thi < 0
            phi_i = pi;
        end
        dk_rho = k0 * sin(abs(thi));
        phase_OB = -dk_rho * rho_refl .* cos(PH_fo - phi_i);

        int_cen = dot_bs .* exp(1j * phase_OB) .* sin(TH_fo);
        P_cen(ii) = abs(sum(int_cen * dth * dph, 'all'))^2;

        int_disp = dot_d .* exp(1j * phase_OB) .* sin(TH_fo);
        P_disp(ii) = abs(sum(int_disp * dth * dph, 'all'))^2;
    end

    P_norm = max(P_cen);
    P_cen_dB = 10 * log10(P_cen / P_norm);
    P_disp_dB = 10 * log10(P_disp / P_norm);

    % Find beam crossing level
    [~, idx_pk_d] = max(P_disp);
    theta_pk = theta_i_ext(idx_pk_d);

    % Find crossing point between central and displaced beams
    diff_beams = P_cen - P_disp;
    crossings = find(diff_beams(1:end-1) .* diff_beams(2:end) < 0);

    cross_level = NaN;
    cross_angle = NaN;
    if ~isempty(crossings)
        % Find the crossing closest to midpoint between peaks
        mid_angle = theta_pk / 2;
        [~, best_cross_idx] = min(abs(theta_i_ext(crossings) - mid_angle));
        ci = crossings(best_cross_idx);
        % Interpolate
        alpha = abs(diff_beams(ci)) / (abs(diff_beams(ci)) + abs(diff_beams(ci+1)));
        cross_angle = theta_i_ext(ci) + alpha * (theta_i_ext(ci+1) - theta_i_ext(ci));
        cross_power = P_cen(ci) + alpha * (P_cen(ci+1) - P_cen(ci));
        cross_level = 10 * log10(cross_power / P_norm);
    end

    eta_disp = (max(P_disp) / P_norm) * eta_ap_rx;
    theta_s = rad2deg(theta_i_ext(idx_pk_d));

    fprintf('\n=== Case %s: dx = %.2f lambda0 ===\n', case_labels{ic}, dx/lambda0);
    fprintf('Pointing: %.3f deg (theory: %.3f deg)\n', theta_s, rad2deg(asin(dx/F)));
    fprintf('eta_ap: %.4f (%.2f%%)\n', eta_disp, eta_disp*100);
    fprintf('Scan loss: %.2f dB\n', 10*log10(max(P_disp)/P_norm));
    fprintf('Beam crossing level: %.2f dB at theta = %.3f deg\n', cross_level, rad2deg(cross_angle));

    subplot(1,3,ic);
    plot(rad2deg(theta_i_ext), P_cen_dB, 'b-', 'LineWidth', 1.5); hold on;
    plot(rad2deg(theta_i_ext), P_disp_dB, '-', 'Color', colors{ic}, 'LineWidth', 1.5);
    if ~isnan(cross_level)
        yline(cross_level, 'k--', sprintf('%.1f dB', cross_level), 'LineWidth', 1);
    end
    xlabel('\theta_i [deg]');
    ylabel('Normalized P_{rx} [dB]');
    title(sprintf('d_x = %s = %.1f\\lambda_0', case_labels{ic}, dx/lambda0));
    legend('Central', 'Displaced', 'Location', 'south');
    grid on; ylim([-30 0]);
    xlim([-2 4]);
    set(gca, 'Color', 'w');
end

sgtitle('Q2.2: Sampling Schemes — Beam Crossing Levels');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q3_sampling_schemes.png'));

%% Summary table
fprintf('\n=== Sampling Summary ===\n');
fprintf('Scheme             | dx/lambda0 | Theor. cross | Obtained | Notes\n');
fprintf('Max-gain (2lam f#) | %.1f       | -10 dB       | computed | beams at 2*HPBW apart\n', dx_cases(1)/lambda0);
fprintf('Field (lam f#)     | %.1f       | -3 dB        | computed | beams at HPBW apart\n', dx_cases(2)/lambda0);
fprintf('Power (0.5lam f#)  | %.1f       | -1.5 dB      | computed | beams at HPBW/2 apart\n', dx_cases(3)/lambda0);
fprintf('\nDiscrepancies arise because theoretical values assume an Airy (uniform) beam.\n');
fprintf('The actual beam is broadened by the feed taper (~-14 dB edge illumination),\n');
fprintf('causing beams to cross at higher levels for the same angular separation.\n');
