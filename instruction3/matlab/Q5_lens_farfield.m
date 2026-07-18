%% Q5: Lens Far Field vs Airy Pattern

clear; close all; clc;
set_plot_defaults();

c = 3e8; freq = 200e9; lambda0 = c / freq;
k0 = 2*pi/lambda0; er = 4; zeta_0 = 120*pi;
D = 9*lambda0; theta_0 = 40*pi/180;
u0 = 0.5; v0 = 0.5; a_ap = D/2;

lens = ellipse_geometry(er, D, lambda0);
e = lens.e; a = lens.a;

% Aperture grid
N = 301;
x_ap = linspace(-a_ap, a_ap, N);
y_ap = linspace(-a_ap, a_ap, N);
dx = x_ap(2)-x_ap(1); dy = y_ap(2)-y_ap(1);
[X_ap, Y_ap] = meshgrid(x_ap, y_ap);
rho_ap = sqrt(X_ap.^2 + Y_ap.^2);
phi_ap = atan2(Y_ap, X_ap);
phi_ap(rho_ap < 1e-15) = 0;

% Map rho -> theta via interpolation
theta_grid = zeros(size(rho_ap));
theta_fine = linspace(0, theta_0, 10001);
r_fine = a*(1-e^2)./(1-e*cos(theta_fine));
rho_fine = r_fine.*sin(theta_fine);
for i = 1:N
    for j = 1:N
        if rho_ap(i,j) <= a_ap && rho_ap(i,j) > 0
            theta_grid(i,j) = interp1(rho_fine, theta_fine, rho_ap(i,j), 'linear', theta_0);
        end
    end
end

r_at_theta0 = a*(1-e^2)/(1-e*cos(theta_0));
rho_max = r_at_theta0 * sin(theta_0);
mask = rho_ap <= rho_max;

r_ell = a*(1-e^2)./(1-e*cos(theta_grid));
r_ell(~mask) = 1;

theta_i = lens_incidence_angle(e, theta_grid);
[tau_perp, tau_par, ~, ~] = fresnel_coeff(er, theta_i);

u_grid = sin(theta_grid).*cos(phi_ap);
v_grid = sin(theta_grid).*sin(phi_ap);
gauss = exp(-((u_grid/u0).^2 + (v_grid/v0).^2));
Ei_theta = sin(phi_ap).*gauss;
Ei_phi   = cos(phi_ap).*gauss;

Jx = -(2/zeta_0) * (tau_par.*Ei_theta.*cos(phi_ap) - tau_perp.*Ei_phi.*sin(phi_ap)) ./ r_ell;
Jy = -(2/zeta_0) * (tau_par.*Ei_theta.*sin(phi_ap) + tau_perp.*Ei_phi.*cos(phi_ap)) ./ r_ell;
Jx = Jx .* mask;
Jy = Jy .* mask;

% Far field computation
theta_null_airy = asin(1.22*lambda0/D);
theta_max_ff = min(15*theta_null_airy, pi/4);
Nth_ff = 1001;
theta_ff = linspace(0, theta_max_ff, Nth_ff);

E_tot_E = zeros(1,Nth_ff);
E_tot_H = zeros(1,Nth_ff);

for ip = 1:2
    if ip == 1, phi_ff = pi/2; else, phi_ff = 0; end
    for it = 1:Nth_ff
        th = theta_ff(it);
        kx_s = k0*sin(th)*cos(phi_ff);
        ky_s = k0*sin(th)*sin(phi_ff);
        kz_s = k0*cos(th);
        sgf = EJ_SGF(1, k0, kx_s, ky_s);

        phase = exp(1j*(kx_s*X_ap + ky_s*Y_ap));
        Jtx = sum(Jx.*phase,'all')*dx*dy;
        Jty = sum(Jy.*phase,'all')*dx*dy;

        Ex = 1j*kz_s*(sgf.Gxx*Jtx + sgf.Gxy*Jty);
        Ey = 1j*kz_s*(sgf.Gyx*Jtx + sgf.Gyy*Jty);
        Ez = 1j*kz_s*(sgf.Gzx*Jtx + sgf.Gzy*Jty);

        E_th = Ex*cos(th)*cos(phi_ff) + Ey*cos(th)*sin(phi_ff) - Ez*sin(th);
        E_ph = -Ex*sin(phi_ff) + Ey*cos(phi_ff);

        if ip==1, E_tot_E(it) = sqrt(abs(E_th)^2+abs(E_ph)^2);
        else,     E_tot_H(it) = sqrt(abs(E_th)^2+abs(E_ph)^2); end
    end
end

E_max = max([max(E_tot_E), max(E_tot_H)]);
E_E_dB = 20*log10(E_tot_E/E_max + 1e-15);
E_H_dB = 20*log10(E_tot_H/E_max + 1e-15);

% Airy pattern
airy = ones(size(theta_ff));
for it = 1:Nth_ff
    arg = k0*a_ap*sin(theta_ff(it));
    if abs(arg) > 1e-10, airy(it) = 2*besselj(1,arg)/arg; end
end
airy_dB = 20*log10(abs(airy)/max(abs(airy)));

fig = figure('Color','w','Position',[100 100 900 550]);
set(fig,'InvertHardcopy','off');
plot(theta_ff*180/pi, E_E_dB, 'b-', 'LineWidth', 2); hold on;
plot(theta_ff*180/pi, E_H_dB, 'r--', 'LineWidth', 2);
plot(theta_ff*180/pi, airy_dB, 'k:', 'LineWidth', 2);
xlabel('\theta [deg]'); ylabel('Normalized |E^{far}| [dB]');
title('Q5: Lens Far Field (D = 9\lambda_0, \theta_0 = 40°, \epsilon_r = 4)');
lg = legend('E-plane (\phi = 90°)', 'H-plane (\phi = 0°)', 'Airy pattern (uniform)', 'Location', 'northeast');
set(lg, 'FontSize', 10);
grid on; ylim([-40 0]); xlim([0 theta_max_ff*180/pi]);
print('-dpng', '-r150', fullfile('..', 'figures', 'Q5_lens_farfield.png'));
fprintf('Q5 done.\n');
