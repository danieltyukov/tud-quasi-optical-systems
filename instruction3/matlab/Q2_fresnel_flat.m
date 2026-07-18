%% Q2: Fresnel Transmission Coefficients - Flat Interface

clear; close all; clc;
set_plot_defaults();

er = 4; n = sqrt(er);

Ntheta = 1001;
theta_i = linspace(0, pi/2, Ntheta);

[~, ~, Gamma_perp, Gamma_par] = fresnel_coeff(er, theta_i);

Pt_Pi_TE = 1 - abs(Gamma_perp).^2;
Pt_Pi_TM = 1 - abs(Gamma_par).^2;

theta_c = asin(1/n);
Pt_Pi_TE(theta_i > theta_c) = 0;
Pt_Pi_TM(theta_i > theta_c) = 0;

theta_B = atan(1/n);
Pt_Pi_normal = 4*n / (1+n)^2;

fprintf('Critical angle: %.2f deg\n', theta_c * 180/pi);
fprintf('Brewster angle: %.2f deg\n', theta_B * 180/pi);
fprintf('P_t/P_i at normal: %.4f\n', Pt_Pi_normal);

fig = figure('Color','w','Position',[100 100 800 500]);
set(fig,'InvertHardcopy','off');
plot(theta_i*180/pi, Pt_Pi_TE, 'b-', 'LineWidth', 2); hold on;
plot(theta_i*180/pi, Pt_Pi_TM, 'r--', 'LineWidth', 2);
xline(theta_c*180/pi, 'k-.', 'LineWidth', 1.2);
xline(theta_B*180/pi, 'm-.', 'LineWidth', 1.2);
plot(0, Pt_Pi_normal, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
text(2, Pt_Pi_normal - 0.06, sprintf('%.3f', Pt_Pi_normal), 'FontSize', 10);
text(theta_c*180/pi + 1, 0.55, sprintf('\\theta_c = %.1f°', theta_c*180/pi), 'FontSize', 10, 'Color', 'k');
text(theta_B*180/pi - 12, 0.45, sprintf('\\theta_B = %.1f°', theta_B*180/pi), 'FontSize', 10, 'Color', 'm');
xlabel('\theta_i [deg]'); ylabel('P_t / P_i');
title('Q2: Power Transmission - Flat Fused Silica-Air Interface (\epsilon_r = 4)');
lg = legend('TE (perpendicular)', 'TM (parallel)', 'Location', 'east');
set(lg, 'FontSize', 10);
grid on; ylim([0 1.1]); xlim([0 50]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q2_fresnel_flat.png'));
fprintf('Q2 done.\n');
