%% Q3: how many elements are needed to achieve HPBW = 2 deg?
%%
%% Compute the windowed broadside-radiated pattern as a function of array
%% size N x N at the centre frequency, locate the half-power beamwidth,
%% and check the result against the closed-form L = 0.886 * lambda / HPBW
%% rule. Also evaluate band-edge effects (HPBW changes with frequency).

clear; close all; clc;
set_plot_defaults();
load('Q1_design.mat', 'h_best','l_best','w','dx','dy','ZL_best','Mmax', ...
                      'f_lo','f_hi','f_c','lam_c');

c0 = 3e8;

%% Closed-form first estimate
% Uniform line array, broadside, half-power half-beamwidth: 0.886*lambda/L
% (radians). For L = N*dx (centre-to-centre count of N elements)
HPBW_target = deg2rad(2);
L_required  = 0.886*lam_c/HPBW_target;     % aperture length at centre freq
Nx_estimate = ceil(L_required/dx);
fprintf('Closed-form estimate: L = %.1f mm  ->  N = %d elements per dimension\n', ...
        L_required*1e3, Nx_estimate);
fprintf('Band-edge estimate :\n');
for f = [f_lo f_c f_hi]
    L_req_f  = 0.886*(c0/f)/HPBW_target;
    fprintf('  f = %.2f GHz : L >= %.1f mm  (N >= %d)\n', f/1e9, L_req_f*1e3, ceil(L_req_f/dx));
end

%% Numerical: windowed broadside pattern for N = {25, 35, 45, 55, 65}
% E-plane (phi_obs = 0). Broadside scan -> kx0 = ky0 = 0.
N_list  = [25, 35, 45, 55, 65];
theta_obs_deg = linspace(-3, 3, 1801);          % zoom near boresight
theta_obs     = deg2rad(theta_obs_deg);

freq_now = f_c;
k0 = 2*pi*freq_now/c0;

% Active impedance and incident scan current at broadside
Z_act = Z_active_GP(0, 0, k0, dx, dy, w, l_best, h_best, Mmax);
i0    = 1/(ZL_best + Z_act);                    % normalised feed current

% Single-element radiation factor (E-plane phi_obs=0)
kx_obs = k0 * sin(theta_obs);
ky_obs = zeros(size(theta_obs));
G_eff  = EJ_SGF_GP(k0, kx_obs, ky_obs, h_best);
Ikx    = basis_long(kx_obs, k0, l_best);
Jky    = basis_trans(ky_obs, w);
EF     = 1j*k0*cos(theta_obs).*G_eff.Gxx.*Ikx.*Jky;       % (proportional to)

HPBW_sim = nan(numel(N_list), 1);
fig = figure('Color','w','Position',[80 80 1100 460]);
set(fig, 'InvertHardcopy', 'off');

subplot(1,2,1); hold on; grid on;
cmap = lines(numel(N_list));
for kn = 1:numel(N_list)
    N = N_list(kn);
    AFx = sin(N*kx_obs*dx/2)./sin(kx_obs*dx/2);
    AFx(abs(kx_obs*dx/2) < 1e-12) = N;
    AFy = N;                                                 % broadside, ky=0 -> Dirichlet kernel = N
    P   = abs(EF .* AFx .* AFy * i0).^2;
    P   = P / max(P);
    plot(theta_obs_deg, 10*log10(P), '-', 'Color', cmap(kn,:), ...
         'LineWidth', 1.4, 'DisplayName', sprintf('N = %d', N));

    % HPBW: full width where 10*log10(P) > -3 dB
    above_half = find(10*log10(P) > -3);
    if isempty(above_half)
        HPBW_sim(kn) = NaN;
    else
        HPBW_sim(kn) = theta_obs_deg(above_half(end)) - theta_obs_deg(above_half(1));
    end
end
yline(-3, 'k:', '-3 dB', 'HandleVisibility', 'off');
xlabel('\theta_{obs} [deg] (E-plane, broadside)');
ylabel('|E|^2 normalised [dB]');
title('Windowed array pattern at f_c = 29.25 GHz');
legend('Location', 'south'); ylim([-20 0.5]); xlim([-3 3]);

subplot(1,2,2); hold on; grid on;
plot(N_list, HPBW_sim, 'b-o', 'LineWidth', 1.7, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
HPBW_cf = 0.886*(c0/freq_now)./(N_list*dx) * 180/pi;
plot(N_list, HPBW_cf, 'r--', 'LineWidth', 1.5);
yline(2, 'k:', 'target = 2\circ', 'HandleVisibility', 'off');
xlabel('N (elements per dimension)');
ylabel('HPBW [deg]');
title('HPBW scaling: simulation vs 0.886\lambda/L');
legend({'simulation','0.886\lambda/(N d_x)'}, 'Location', 'best');

sgtitle({'Q3: half-power beamwidth of the windowed array, f_c = 29.25 GHz', ...
        sprintf('h = %.2f mm, l = %.2f mm, d_x = d_y = %.2f mm, Z_L = %d \\Omega', ...
                h_best*1e3, l_best*1e3, dx*1e3, ZL_best)}, ...
       'FontSize', 11, 'FontWeight', 'bold');
print('-dpng','-r150', fullfile('..','figures','Q3_HPBW.png'));

% Smallest N that gives HPBW <= 2 deg
N_required_sim = N_list(find(HPBW_sim <= 2, 1, 'first'));
if isempty(N_required_sim), N_required_sim = NaN; end

fprintf('\nNumerical HPBW @ f_c:\n');
for kn = 1:numel(N_list)
    fprintf('  N = %2d : HPBW = %.3f deg\n', N_list(kn), HPBW_sim(kn));
end
fprintf('Smallest N (in the swept set) giving HPBW <= 2 deg: %d\n', N_required_sim);
fprintf('Total elements N x N = %d\n', N_required_sim^2);

%% Pattern at the band-edges with the chosen N (worst case f_lo)
N_chosen = 55;
freqs_now = [f_lo, f_c, f_hi];
HPBW_band = nan(size(freqs_now));
fig2 = figure('Color','w','Position',[80 80 700 420]);
set(fig2, 'InvertHardcopy', 'off');
hold on; grid on;
ec = lines(numel(freqs_now));
for kf = 1:numel(freqs_now)
    k0 = 2*pi*freqs_now(kf)/c0;
    Z_act = Z_active_GP(0, 0, k0, dx, dy, w, l_best, h_best, Mmax);
    i0    = 1/(ZL_best + Z_act);
    kx_obs = k0*sin(theta_obs);
    ky_obs = zeros(size(theta_obs));
    G_eff  = EJ_SGF_GP(k0, kx_obs, ky_obs, h_best);
    Ikx    = basis_long(kx_obs, k0, l_best);
    Jky    = basis_trans(ky_obs, w);
    EF     = 1j*k0*cos(theta_obs).*G_eff.Gxx.*Ikx.*Jky;
    AFx    = sin(N_chosen*kx_obs*dx/2)./sin(kx_obs*dx/2);
    AFx(abs(kx_obs*dx/2) < 1e-12) = N_chosen;
    AFy    = N_chosen;
    P      = abs(EF .* AFx .* AFy * i0).^2; P = P/max(P);
    plot(theta_obs_deg, 10*log10(P), '-', 'Color', ec(kf,:), ...
         'LineWidth', 1.5, 'DisplayName', sprintf('%.2f GHz', freqs_now(kf)/1e9));
    above_half = find(10*log10(P) > -3);
    HPBW_band(kf) = theta_obs_deg(above_half(end)) - theta_obs_deg(above_half(1));
end
yline(-3, 'k:', '-3 dB', 'HandleVisibility', 'off');
xlabel('\theta_{obs} [deg]'); ylabel('|E|^2 normalised [dB]');
title(sprintf('N = %d elements: pattern across the band', N_chosen));
legend('Location', 'south'); ylim([-20 0.5]); xlim([-3 3]);
print('-dpng','-r150', fullfile('..','figures','Q3_HPBW_band.png'));

fprintf('\nBand-edge HPBW with N = %d:\n', N_chosen);
for kf = 1:numel(freqs_now)
    fprintf('  f = %.2f GHz: HPBW = %.3f deg\n', freqs_now(kf)/1e9, HPBW_band(kf));
end

save('Q3_HPBW.mat', 'N_list','HPBW_sim','HPBW_cf','N_required_sim','N_chosen', ...
     'freqs_now','HPBW_band','theta_obs_deg');
fprintf('\nFigures saved.\n');
