%% Q3: Fresnel Transmission Coefficients - Elliptical Lens

clear; close all; clc;
set_plot_defaults();

c = 3e8; freq = 200e9; lambda0 = c / freq;
er = 4; D = 9 * lambda0;

lens = ellipse_geometry(er, D, lambda0);
fprintf('Lens: e=%.3f, a=%.4f mm, b=%.4f mm, theta_max=%.1f deg\n', ...
        lens.e, lens.a*1e3, lens.b*1e3, lens.theta_max*180/pi);

Ntheta = 1001;
theta = linspace(0, lens.theta_max, Ntheta);
theta_i_lens = lens_incidence_angle(lens.e, theta);
[~, ~, Gamma_perp, Gamma_par] = fresnel_coeff(er, theta_i_lens);
Pt_Pi_TE = 1 - abs(Gamma_perp).^2;
Pt_Pi_TM = 1 - abs(Gamma_par).^2;

% Flat interface for comparison
theta_c_flat = asin(1/sqrt(er));
theta_i_flat = linspace(0, pi/2, Ntheta);
[~, ~, Gp_flat, Gpar_flat] = fresnel_coeff(er, theta_i_flat);
Pt_TE_flat = 1 - abs(Gp_flat).^2;
Pt_TM_flat = 1 - abs(Gpar_flat).^2;
Pt_TE_flat(theta_i_flat > theta_c_flat) = 0;
Pt_TM_flat(theta_i_flat > theta_c_flat) = 0;

fprintf('theta_i at edge: %.2f deg\n', theta_i_lens(end)*180/pi);

fig = figure('Color','w','Position',[100 100 900 550]);
set(fig,'InvertHardcopy','off');
plot(theta*180/pi, Pt_Pi_TE, 'b-', 'LineWidth', 2.5); hold on;
plot(theta*180/pi, Pt_Pi_TM, 'r-', 'LineWidth', 2.5);
plot(theta_i_flat*180/pi, Pt_TE_flat, 'b--', 'LineWidth', 1.5);
plot(theta_i_flat*180/pi, Pt_TM_flat, 'r--', 'LineWidth', 1.5);
xline(theta_c_flat*180/pi, 'k-.', 'LineWidth', 1.2);
xline(lens.theta_max*180/pi, 'g-.', 'LineWidth', 1.2);
text(theta_c_flat*180/pi + 1, 0.15, sprintf('\\theta_c (flat) = %.0f°', theta_c_flat*180/pi), 'FontSize', 10);
text(lens.theta_max*180/pi + 1, 0.15, sprintf('\\theta_{max} = %.0f°', lens.theta_max*180/pi), 'FontSize', 10, 'Color', [0 0.5 0]);
xlabel('\theta [deg]'); ylabel('P_t / P_i');
title('Q3: Power Transmission - Elliptical Lens vs Flat Interface (\epsilon_r = 4, D = 9\lambda_0)');
lg = legend('TE - Lens', 'TM - Lens', 'TE - Flat', 'TM - Flat', 'Location', 'southwest');
set(lg, 'FontSize', 10);
grid on; ylim([0 1.1]); xlim([0 75]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q3_fresnel_lens.png'));
fprintf('Q3 done.\n');
