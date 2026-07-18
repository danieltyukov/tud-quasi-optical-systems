%% Q1: Feed Far Field in Infinite Dielectric

clear; close all; clc;
set_plot_defaults();

c = 3e8; freq = 200e9; lambda0 = c / freq;
er = 4; u0 = 0.5; v0 = 0.5;

Ntheta = 901;
theta = linspace(0, pi/2, Ntheta);

% E-plane (phi=90): v = sin(theta), u = 0
gauss_E = exp(-(sin(theta)/v0).^2);
% H-plane (phi=0): u = sin(theta), v = 0
gauss_H = exp(-(sin(theta)/u0).^2);

E_E_dB = 20 * log10(gauss_E + 1e-15);
E_H_dB = 20 * log10(gauss_H + 1e-15);

theta_10dB = asin(min(u0*sqrt(log(10)), 1));
theta_3dB  = asin(u0*sqrt(log(2)));
fprintf('Feed -3 dB half-angle:  %.2f deg\n', theta_3dB * 180/pi);
fprintf('Feed -10 dB half-angle: %.2f deg\n', theta_10dB * 180/pi);

fig = figure('Color','w','Position',[100 100 800 500]);
set(fig,'InvertHardcopy','off');
plot(theta*180/pi, E_E_dB, 'b-', 'LineWidth', 2); hold on;
plot(theta*180/pi, E_H_dB, 'r--', 'LineWidth', 2);
xline(theta_10dB*180/pi, 'k:', 'LineWidth', 1);
text(theta_10dB*180/pi + 1, -11, sprintf('-10 dB at %.1f°', theta_10dB*180/pi), 'FontSize', 10);
xlabel('\theta [deg]'); ylabel('Normalized |E^{far}| [dB]');
title('Q1: Feed Far Field in Fused Silica (\epsilon_r = 4, u_0 = v_0 = 0.5)');
lg = legend('E-plane (\phi = 90°)', 'H-plane (\phi = 0°)', 'Location', 'southwest');
set(lg, 'FontSize', 10);
grid on; ylim([-40 0]); xlim([0 90]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q1_feed_farfield.png'));
fprintf('Q1 done.\n');
