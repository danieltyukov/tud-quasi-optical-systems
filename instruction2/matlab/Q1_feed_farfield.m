%% Q1: Feed Far Field Patterns (1 point)

clear; close all; clc;
set_plot_defaults();

%% Constants
c = 3e8;
freq = 500e9;
lambda0 = c / freq;
k0 = 2 * pi / lambda0;
D_f = 4 * lambda0;
a_feed = D_f / 2;

%% Compute far field
Ntheta = 901;
theta = linspace(0, pi/2, Ntheta);

% E-plane: phi = pi/2
phi_E = (pi/2) * ones(size(theta));
[E_th_E, E_ph_E] = feed_farfield(k0, a_feed, theta, phi_E);
E_tot_E = sqrt(abs(E_th_E).^2 + abs(E_ph_E).^2);

% H-plane: phi = 0
phi_H = zeros(size(theta));
[E_th_H, E_ph_H] = feed_farfield(k0, a_feed, theta, phi_H);
E_tot_H = sqrt(abs(E_th_H).^2 + abs(E_ph_H).^2);

% Normalize
E_max = max([max(E_tot_E), max(E_tot_H)]);
E_E_dB = 20 * log10(E_tot_E / E_max);
E_H_dB = 20 * log10(E_tot_H / E_max);

theta_null = asin(1.22 * lambda0 / D_f);
fprintf('Theoretical first null at: %.2f deg\n', theta_null * 180/pi);

%% Plot
fig = figure('Color', 'w', 'Position', [100 100 800 500]);
set(fig, 'InvertHardcopy', 'off');
plot(theta * 180/pi, E_E_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta * 180/pi, E_H_dB, 'r--', 'LineWidth', 1.5);
xline(theta_null * 180/pi, 'k:', sprintf('First null = %.1f°', theta_null*180/pi), ...
      'LineWidth', 1, 'LabelOrientation', 'horizontal');
xlabel('\theta [deg]');
ylabel('Normalized |E^{far}| [dB]');
title('Q1: Feed Far Field - Circular Aperture (D_f = 4\lambda, 500 GHz, y-polarized)');
legend('E-plane (\phi = 90°)', 'H-plane (\phi = 0°)', 'Location', 'southwest');
grid on;
ylim([-40 0]);
xlim([0 90]);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_feed_farfield.png'));
