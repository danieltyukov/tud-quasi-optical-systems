%% Q1.1 + Q1.2: GO field and feed far field on the FO sphere
clear; close all; clc;
set_plot_defaults();

%% Parameters
c = 3e8;
freq = 220e9;
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
zeta0 = 120 * pi;
E0_PW = 1; % V/m

D_f = 3.6 * lambda0;
a_feed = D_f / 2;
D_r = 130 * lambda0;
f_num = 1.8;
F = f_num * D_r;
theta0 = 2 * atan(D_r / (4 * F));

fprintf('lambda0 = %.4f mm\n', lambda0*1e3);
fprintf('D_f = %.4f mm (%.1f lambda)\n', D_f*1e3, D_f/lambda0);
fprintf('D_r = %.4f mm (%.1f lambda)\n', D_r*1e3, D_r/lambda0);
fprintf('F = %.4f mm\n', F*1e3);
fprintf('theta0 = %.2f deg\n', rad2deg(theta0));

%% Compute fields on FO sphere
Nth = 501;
theta = linspace(0, theta0, Nth);

% --- phi = 0 plane (E-plane for x-pol) ---
phi_E = zeros(size(theta));
[Eth_feed_E, Eph_feed_E] = feed_farfield_x(k0, a_feed, theta, phi_E);
Etot_feed_E = sqrt(abs(Eth_feed_E).^2 + abs(Eph_feed_E).^2);

% GO field: E_theta = -2*cos(phi)/(1+cos(theta)), E_phi = 2*sin(phi)/(1+cos(theta))
S = 2 ./ (1 + cos(theta));
Eth_GO_E = -S .* cos(phi_E) * E0_PW;  % = -S
Eph_GO_E = S .* sin(phi_E) * E0_PW;   % = 0
Etot_GO_E = sqrt(abs(Eth_GO_E).^2 + abs(Eph_GO_E).^2);

% --- phi = 90 plane (H-plane for x-pol) ---
phi_H = (pi/2) * ones(size(theta));
[Eth_feed_H, Eph_feed_H] = feed_farfield_x(k0, a_feed, theta, phi_H);
Etot_feed_H = sqrt(abs(Eth_feed_H).^2 + abs(Eph_feed_H).^2);

Eth_GO_H = -S .* cos(phi_H) * E0_PW;  % = 0
Eph_GO_H = S .* sin(phi_H) * E0_PW;   % = S
Etot_GO_H = sqrt(abs(Eth_GO_H).^2 + abs(Eph_GO_H).^2);

%% Normalize feed far field to match GO at theta=0
feed_norm = Etot_feed_E(1);
GO_norm = Etot_GO_E(1);
scale = GO_norm / feed_norm;

Etot_feed_E_scaled = Etot_feed_E * scale;
Etot_feed_H_scaled = Etot_feed_H * scale;

%% Plot
fig = figure('Color', 'w', 'Position', [100 100 1000 450]);
set(fig, 'InvertHardcopy', 'off');

subplot(1,2,1);
plot(rad2deg(theta), Etot_GO_E, 'b-', 'LineWidth', 2); hold on;
plot(rad2deg(theta), Etot_feed_E_scaled, 'r--', 'LineWidth', 1.5);
xlabel('\theta [deg]');
ylabel('|E| [V/m]');
title('\phi = 0° (E-plane)');
legend('GO field', 'Feed far field (scaled)', 'Location', 'best');
grid on;
xlim([0 rad2deg(theta0)]);
set(gca, 'Color', 'w');

subplot(1,2,2);
plot(rad2deg(theta), Etot_GO_H, 'b-', 'LineWidth', 2); hold on;
plot(rad2deg(theta), Etot_feed_H_scaled, 'r--', 'LineWidth', 1.5);
xlabel('\theta [deg]');
ylabel('|E| [V/m]');
title('\phi = 90° (H-plane)');
legend('GO field', 'Feed far field (scaled)', 'Location', 'best');
grid on;
xlim([0 rad2deg(theta0)]);
set(gca, 'Color', 'w');

sgtitle('Q1: GO Field vs Feed Far Field on FO Sphere (220 GHz, f_# = 1.8)');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_GO_and_feed.png'));

%% Also plot in dB (unnormalized comparison)
fig2 = figure('Color', 'w', 'Position', [100 100 1000 450]);
set(fig2, 'InvertHardcopy', 'off');

subplot(1,2,1);
plot(rad2deg(theta), 20*log10(Etot_GO_E/max(Etot_GO_E)), 'b-', 'LineWidth', 2); hold on;
plot(rad2deg(theta), 20*log10(Etot_feed_E_scaled/max(Etot_feed_E_scaled)), 'r--', 'LineWidth', 1.5);
xlabel('\theta [deg]');
ylabel('Normalized |E| [dB]');
title('\phi = 0° (E-plane)');
legend('GO field', 'Feed far field', 'Location', 'southwest');
grid on; ylim([-30 0]);
xlim([0 rad2deg(theta0)]);
set(gca, 'Color', 'w');

subplot(1,2,2);
plot(rad2deg(theta), 20*log10(Etot_GO_H/max(Etot_GO_H)), 'b-', 'LineWidth', 2); hold on;
plot(rad2deg(theta), 20*log10(Etot_feed_H_scaled/max(Etot_feed_H_scaled)), 'r--', 'LineWidth', 1.5);
xlabel('\theta [deg]');
ylabel('Normalized |E| [dB]');
title('\phi = 90° (H-plane)');
legend('GO field', 'Feed far field', 'Location', 'southwest');
grid on; ylim([-30 0]);
xlim([0 rad2deg(theta0)]);
set(gca, 'Color', 'w');

sgtitle('Q1: GO vs Feed Far Field [dB] on FO Sphere');
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_GO_and_feed_dB.png'));

fprintf('\nEdge taper (feed at theta0):\n');
fprintf('  E-plane: %.1f dB\n', 20*log10(Etot_feed_E_scaled(end)/Etot_feed_E_scaled(1)));
fprintf('  H-plane: %.1f dB\n', 20*log10(Etot_feed_H_scaled(end)/Etot_feed_H_scaled(1)));
