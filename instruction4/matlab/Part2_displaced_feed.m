%% Assignment 4 - Part 2: Feed displaced in the focal plane
%  FO for FPA - Displaced Feed Analysis
%
%  The circular feed is displaced by d_x in the focal plane.
%  In Q2.2, both d_x AND D_f change for each sampling scenario.

clear; close all; clc;
set_plot_defaults();

%% Load Part 1 workspace (or recompute if not available)
if exist(fullfile('..', 'matlab', 'Part1_workspace.mat'), 'file')
    load(fullfile('..', 'matlab', 'Part1_workspace.mat'));
    fprintf('Loaded Part 1 workspace.\n');
else
    error('Run Part1_FO_FPA.m first to generate the workspace.');
end

fprintf('=== Assignment 4 - Part 2: Displaced Feed ===\n');
fprintf('Reflector: D_r = %.0f lambda0, f# = %.1f, F = %.4f m\n', D_r/lambda0, f_num, F);

%% ========== Q2.1: d_x = D_f = 3.6*lambda0 ==========
fprintf('\n--- Q2.1: d_x = D_f = 3.6*lambda0 ---\n');

d_x = D_f;  % = 3.6*lambda0
fprintf('d_x = %.2f lambda0 = %.4f mm\n', d_x/lambda0, d_x*1e3);

% Pointing direction: sin(theta_s) = d_x / F
theta_s = asin(d_x / F);
fprintf('Pointing direction: theta_s = %.5f deg\n', theta_s*180/pi);

% Angular grid (loaded from Part 1)
% TH, PH are Nph x Nth meshgrids

% Phase shift for displaced feed (in Tx far-field domain)
% Displacing feed by d_x along x introduces phase: exp(j*k0*sin(theta)*cos(phi)*d_x)
phase_shift = exp(1j * k0 * sin(TH) .* cos(PH) * d_x);
Va_th_disp = Va_th .* phase_shift;
Va_ph_disp = Va_ph .* phase_shift;

% Compute Rx pattern for displaced feed over theta_i range
% Use a finer theta_i range for better pattern resolution
N_thi_q21 = 91;
theta_i_q21 = linspace(0, 3*pi/180, N_thi_q21);

% Centered beam pattern (phi_i=0 cut)
P_rx_center = zeros(N_thi_q21, 1);
P_rx_disp_q21 = zeros(N_thi_q21, 1);

fprintf('Computing centered and displaced Rx patterns...\n');
for i_thi = 1:N_thi_q21
    thi = theta_i_q21(i_thi);

    % GO field for this incidence angle
    [eGO_th_off, eGO_ph_off] = GO_field_off_broadside(E0_PW, k0, F, ...
                                    TH, PH, theta_0, thi, 0);

    % Centered feed
    Voc_c = compute_Voc(Va_th, Va_ph, eGO_th_off, eGO_ph_off, TH, dth, dph, zeta);
    P_rx_center(i_thi) = abs(Voc_c)^2 / (16 * P_rad_feed);

    % Displaced feed
    Voc_d = compute_Voc(Va_th_disp, Va_ph_disp, eGO_th_off, eGO_ph_off, TH, dth, dph, zeta);
    P_rx_disp_q21(i_thi) = abs(Voc_d)^2 / (16 * P_rad_feed);
end

% Aperture efficiency
P_inc = (abs(E0_PW)^2 / (2*zeta)) * pi * a_refl^2;
eta_ap_center = max(P_rx_center) / P_inc;
eta_ap_disp = max(P_rx_disp_q21) / P_inc;

fprintf('Aperture efficiency (centered)  = %.4f (%.2f%%)\n', eta_ap_center, eta_ap_center*100);
fprintf('Aperture efficiency (displaced) = %.4f (%.2f%%)\n', eta_ap_disp, eta_ap_disp*100);

% Normalize to centered beam max
P_max_ref = max(P_rx_center);
P_center_dB = 10*log10(P_rx_center / P_max_ref + 1e-30);
P_disp_dB = 10*log10(P_rx_disp_q21 / P_max_ref + 1e-30);

fig1 = figure('Color', 'w', 'Position', [100 100 900 500]);
plot(theta_i_q21*180/pi, P_center_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_i_q21*180/pi, P_disp_dB, 'r--', 'LineWidth', 1.5);
xline(theta_s*180/pi, 'k:', sprintf('\\theta_s = %.3f°', theta_s*180/pi), 'LineWidth', 1);
xlabel('\theta_i [deg]');
ylabel('Normalized Rx power [dB]');
title(sprintf('Q2.1: Centered vs Displaced (d_x = %.1f\\lambda_0, D_f = %.1f\\lambda_0)', ...
      d_x/lambda0, D_f/lambda0));
legend('Centered feed', sprintf('Displaced d_x = %.1f\\lambda_0', d_x/lambda0), ...
       'Location', 'best');
grid on; ylim([-40 0]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_1_centered_vs_displaced.png'));

fprintf('\nDifferences between centered and displaced beams:\n');
fprintf('  - Displaced beam points at theta_s = %.5f deg (not broadside)\n', theta_s*180/pi);
fprintf('  - The two beams should have similar shapes for large f/#\n');
fprintf('  - Any asymmetry in the displaced beam is due to coma aberration\n');
fprintf('  - For f/# = %.1f, coma should be small\n', f_num);

%% ========== Q2.2: Three sampling scenarios ==========
fprintf('\n\n--- Q2.2: Three sampling scenarios ---\n');

% Assignment says d_x = D_f for each case. This means both the feed
% diameter and the displacement change together:
% a) d_x = D_f = 2*lambda0*f#     = 3.6*lambda0  (max gain sampling)
% b) d_x = D_f = lambda0*f#       = 1.8*lambda0  (field sampling)
% c) d_x = D_f = 0.5*lambda0*f#   = 0.9*lambda0  (power sampling)

d_cases = [2*lambda0*f_num, lambda0*f_num, 0.5*lambda0*f_num];
Df_cases = d_cases;  % D_f = d_x in each case
case_names = {'a) Max gain (2\lambda_0 f_#)', ...
              'b) Field (\lambda_0 f_#)', ...
              'c) Power (0.5\lambda_0 f_#)'};
case_short = {'max gain', 'field', 'power'};

% Theoretical beam crossing levels from lecture notes:
% Max gain sampling: beams cross at ~-10 dB
% Field sampling:    beams cross at ~-3 dB
% Power sampling:    beams cross at ~-1.5 dB
theoretical_crossing_dB = [-10, -3, -1.5];

% Angular scanning separation
% theta_HPBW ~ lambda0 / D_r
theta_HPBW = lambda0 / D_r;  % in radians
fprintf('theta_HPBW ~ lambda0/D_r = %.5f deg\n', theta_HPBW*180/pi);

% Compute patterns for each case
N_thi_ext = 151;
theta_i_ext = linspace(0, 3*pi/180, N_thi_ext);

fig2 = figure('Color', 'w', 'Position', [100 100 1100 700]);
colors = {'r', [0 0.6 0], [0.8 0 0.8]};
line_styles = {'-', '--', '-.'};

for ic = 1:3
    d_x_case = d_cases(ic);
    Df_case = Df_cases(ic);
    a_feed_case = Df_case / 2;

    fprintf('\nCase %c: d_x = D_f = %.3f lambda0 = %.2f*lambda0*f#\n', ...
            char('a'+ic-1), d_x_case/lambda0, d_x_case/(lambda0*f_num));

    % Recompute feed far field with new feed diameter
    [Va_th_case, Va_ph_case] = feed_farfield_x(k0, a_feed_case, TH, PH);

    % Feed radiated power for this feed size
    U_feed_case = abs(Va_th_case).^2 + abs(Va_ph_case).^2;
    P_rad_feed_case = (1/zeta) * sum(U_feed_case(:) .* sin(TH(:))) * dth * dph;

    % Displaced feed: add phase shift
    phase_shift_case = exp(1j * k0 * sin(TH) .* cos(PH) * d_x_case);
    Va_th_d = Va_th_case .* phase_shift_case;
    Va_ph_d = Va_ph_case .* phase_shift_case;

    % Pointing direction
    theta_s_case = asin(d_x_case / F);
    fprintf('  Pointing: theta_s = %.5f deg\n', theta_s_case*180/pi);
    fprintf('  theta_s / theta_HPBW = %.2f\n', theta_s_case / theta_HPBW);

    % Compute Rx patterns: centered + displaced (phi_i = 0 cut)
    P_rx_c = zeros(N_thi_ext, 1);
    P_rx_d = zeros(N_thi_ext, 1);

    for i_thi = 1:N_thi_ext
        thi = theta_i_ext(i_thi);
        [eGO_th_off, eGO_ph_off] = GO_field_off_broadside(E0_PW, k0, F, ...
                                        TH, PH, theta_0, thi, 0);

        % Centered feed (this feed size)
        Voc_c = compute_Voc(Va_th_case, Va_ph_case, eGO_th_off, eGO_ph_off, ...
                            TH, dth, dph, zeta);
        P_rx_c(i_thi) = abs(Voc_c)^2 / (16 * P_rad_feed_case);

        % Displaced feed
        Voc_d = compute_Voc(Va_th_d, Va_ph_d, eGO_th_off, eGO_ph_off, ...
                            TH, dth, dph, zeta);
        P_rx_d(i_thi) = abs(Voc_d)^2 / (16 * P_rad_feed_case);
    end

    % Aperture efficiency
    eta_ap_c_case = max(P_rx_c) / P_inc;
    eta_ap_d_case = max(P_rx_d) / P_inc;
    fprintf('  eta_ap (centered, D_f=%.1f lam) = %.4f\n', Df_case/lambda0, eta_ap_c_case);
    fprintf('  eta_ap (displaced) = %.4f\n', eta_ap_d_case);

    % Normalize to centered beam max for this feed size
    P_max_c = max(P_rx_c);
    P_c_dB = 10*log10(P_rx_c / P_max_c + 1e-30);
    P_d_dB = 10*log10(P_rx_d / P_max_c + 1e-30);

    % Find beam crossing level
    % The crossing happens where P_rx_center = P_rx_displaced
    diff_beams = P_c_dB - P_d_dB;
    sign_changes = find(diff(sign(diff_beams)) ~= 0);
    crossing_level_actual = NaN;
    crossing_theta = NaN;
    if ~isempty(sign_changes)
        idx_c = sign_changes(1);
        x1 = theta_i_ext(idx_c)*180/pi; x2 = theta_i_ext(idx_c+1)*180/pi;
        y1 = diff_beams(idx_c); y2 = diff_beams(idx_c+1);
        crossing_theta = x1 - y1*(x2-x1)/(y2-y1);
        crossing_level_actual = interp1(theta_i_ext*180/pi, P_c_dB, crossing_theta);
    end
    fprintf('  Beam crossing: theta = %.4f deg, level = %.2f dB (theory: %.1f dB)\n', ...
            crossing_theta, crossing_level_actual, theoretical_crossing_dB(ic));

    if abs(crossing_level_actual - theoretical_crossing_dB(ic)) > 1.5
        fprintf('  >> DISCREPANCY: Actual crossing differs from theory.\n');
        fprintf('     This is because the feed pattern is not a Dirac delta;\n');
        fprintf('     the finite feed size affects the effective beamwidth\n');
        fprintf('     and thus the crossing level.\n');
    end

    % Plot
    subplot(1,3,ic);
    plot(theta_i_ext*180/pi, P_c_dB, 'b-', 'LineWidth', 1.5); hold on;
    plot(theta_i_ext*180/pi, P_d_dB, '--', 'Color', colors{ic}, 'LineWidth', 1.5);

    % Mark crossing point
    if ~isnan(crossing_theta)
        plot(crossing_theta, crossing_level_actual, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
        text(crossing_theta + 0.05, crossing_level_actual - 2, ...
             sprintf('%.1f dB', crossing_level_actual), 'FontSize', 9);
    end

    yline(theoretical_crossing_dB(ic), 'k:', sprintf('Theory: %.1f dB', ...
          theoretical_crossing_dB(ic)), 'LineWidth', 0.8);

    xlabel('\theta_i [deg]');
    ylabel('Normalized power [dB]');
    title(sprintf('Case %c: D_f = %.1f\\lambda_0', char('a'+ic-1), Df_case/lambda0));
    legend('Centered', 'Displaced', 'Location', 'best');
    grid on; ylim([-25 0]); xlim([0 2.5]);
end

sgtitle('Q2.2: Beam Crossing Levels for Three FPA Sampling Scenarios');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_2_sampling_scenarios.png'));

%% Summary
fprintf('\n\n============ SUMMARY TABLE ============\n');
fprintf('%-20s  %8s  %8s  %10s  %10s  %12s\n', ...
        'Sampling', 'd_x/lam', 'D_f/lam', 'theta_s', 'eta_ap', 'Crossing dB');
fprintf('%s\n', repmat('-', 1, 75));

for ic = 1:3
    d_x_case = d_cases(ic);
    Df_case = Df_cases(ic);
    a_feed_case = Df_case / 2;
    theta_s_case = asin(d_x_case / F);

    % Recompute efficiency for displaced feed at pointing direction
    [Va_th_case, Va_ph_case] = feed_farfield_x(k0, a_feed_case, TH, PH);
    U_feed_case = abs(Va_th_case).^2 + abs(Va_ph_case).^2;
    P_rad_case = (1/zeta) * sum(U_feed_case(:) .* sin(TH(:))) * dth * dph;

    phase_shift_case = exp(1j * k0 * sin(TH) .* cos(PH) * d_x_case);
    Va_th_d = Va_th_case .* phase_shift_case;
    Va_ph_d = Va_ph_case .* phase_shift_case;

    [eGO_pt_th, eGO_pt_ph] = GO_field_off_broadside(E0_PW, k0, F, ...
                                    TH, PH, theta_0, theta_s_case, 0);
    Voc_pt = compute_Voc(Va_th_d, Va_ph_d, eGO_pt_th, eGO_pt_ph, TH, dth, dph, zeta);
    P_rx_pt = abs(Voc_pt)^2 / (16 * P_rad_case);
    eta_case = P_rx_pt / P_inc;

    % Centered beam eta_ap for this feed size
    [eGO_bs_th, eGO_bs_ph] = GO_field_parabolic(E0_PW, TH, PH, theta_0);
    Voc_bs = compute_Voc(Va_th_case, Va_ph_case, eGO_bs_th, eGO_bs_ph, TH, dth, dph, zeta);
    eta_center = abs(Voc_bs)^2 / (16 * P_rad_case) / P_inc;

    fprintf('%-20s  %8.2f  %8.2f  %10.4f°  %10.4f  %10.1f (th)\n', ...
            case_short{ic}, d_x_case/lambda0, Df_case/lambda0, ...
            theta_s_case*180/pi, eta_center, theoretical_crossing_dB(ic));
end
fprintf('%s\n', repmat('=', 1, 75));

fprintf('\nKey observations:\n');
fprintf('  - Max gain sampling (2*lam*f#): largest feed, highest efficiency, sparsest coverage\n');
fprintf('  - Field sampling (lam*f#):      beams cross at -3 dB, Nyquist-like coverage\n');
fprintf('  - Power sampling (0.5*lam*f#):  smallest feed, lowest efficiency, densest coverage\n');
fprintf('  - The efficiency vs FoV sampling trade-off is independent of f/#\n');
fprintf('  - Discrepancies from theoretical crossing levels arise because the\n');
fprintf('    feed has a non-ideal (non-delta) pattern that modifies the effective\n');
fprintf('    system beamwidth compared to the pure Airy pattern assumption.\n');

fprintf('\nPart 2 complete.\n');
