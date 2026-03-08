%% Q6: Directivity and Gain Estimation

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
rho_max = r_at_theta0*sin(theta_0);
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

Jx = -(2/zeta_0)*(tau_par.*Ei_theta.*cos(phi_ap) - tau_perp.*Ei_phi.*sin(phi_ap))./r_ell;
Jy = -(2/zeta_0)*(tau_par.*Ei_theta.*sin(phi_ap) + tau_perp.*Ei_phi.*cos(phi_ap))./r_ell;
Jx = Jx.*mask; Jy = Jy.*mask;

% Far field on 2D grid for directivity
Nth_ff = 301; Nph_ff = 181;
theta_ff = linspace(0, pi/2 - 1e-4, Nth_ff);
phi_ff = linspace(0, 2*pi, Nph_ff);
dth_ff = theta_ff(2)-theta_ff(1);
dph_ff = phi_ff(2)-phi_ff(1);
U = zeros(Nph_ff, Nth_ff);

fprintf('Computing far-field pattern for directivity...\n');
for ip = 1:Nph_ff
    for it = 1:Nth_ff
        th = theta_ff(it); ph = phi_ff(ip);
        kx_s = k0*sin(th)*cos(ph);
        ky_s = k0*sin(th)*sin(ph);
        kz_s = k0*cos(th);
        sgf = EJ_SGF(1, k0, kx_s, ky_s);

        phase = exp(1j*(kx_s*X_ap + ky_s*Y_ap));
        Jtx = sum(Jx.*phase,'all')*dx*dy;
        Jty = sum(Jy.*phase,'all')*dx*dy;

        vals = 1j*kz_s*[sgf.Gxx sgf.Gxy sgf.Gyx sgf.Gyy sgf.Gzx sgf.Gzy];
        vals(isnan(vals)|isinf(vals)) = 0;

        Ex = vals(1)*Jtx + vals(2)*Jty;
        Ey = vals(3)*Jtx + vals(4)*Jty;
        Ez = vals(5)*Jtx + vals(6)*Jty;

        E_th = Ex*cos(th)*cos(ph) + Ey*cos(th)*sin(ph) - Ez*sin(th);
        E_ph = -Ex*sin(ph) + Ey*cos(ph);
        U(ip,it) = abs(E_th)^2 + abs(E_ph)^2;
    end
end

U_max = max(U(:));
P_rad = sum(U .* sin(theta_ff) * dth_ff * dph_ff, 'all');
D_val = 4*pi*U_max / P_rad;
D_dBi = 10*log10(D_val);
D_max_val = (pi*D/lambda0)^2;
D_max_dBi = 10*log10(D_max_val);

% Reflection efficiency
Nth_f = 901; Nph_f = 361;
theta_feed = linspace(0, pi/2, Nth_f);
phi_feed = linspace(0, 2*pi, Nph_f);
dth_f = theta_feed(2)-theta_feed(1);
dph_f = phi_feed(2)-phi_feed(1);
[TH_f, PH_f] = meshgrid(theta_feed, phi_feed);

u_f = sin(TH_f).*cos(PH_f); v_f = sin(TH_f).*sin(PH_f);
U_feed = exp(-2*((u_f/u0).^2 + (v_f/v0).^2));
P_feed_total = sum(U_feed .* sin(TH_f) * dth_f * dph_f, 'all');

idx_th0 = find(theta_feed <= theta_0, 1, 'last');
P_intercepted = sum(U_feed(:,1:idx_th0) .* sin(TH_f(:,1:idx_th0)) * dth_f * dph_f, 'all');

theta_i_f = lens_incidence_angle(e, TH_f(:,1:idx_th0));
[~, ~, Gp_f, Gpar_f] = fresnel_coeff(er, theta_i_f);
gauss_f = exp(-((u_f(:,1:idx_th0)/u0).^2 + (v_f(:,1:idx_th0)/v0).^2));
Ei_par_f = sin(PH_f(:,1:idx_th0)).*gauss_f;
Ei_perp_f = cos(PH_f(:,1:idx_th0)).*gauss_f;
P_ref = abs(Gpar_f).^2.*abs(Ei_par_f).^2 + abs(Gp_f).^2.*abs(Ei_perp_f).^2;
P_inc = abs(Ei_par_f).^2 + abs(Ei_perp_f).^2;
P_transmitted = sum((P_inc - P_ref) .* sin(TH_f(:,1:idx_th0)) * dth_f * dph_f, 'all');

eta_spill = P_intercepted / P_feed_total;
eta_refl = P_transmitted / P_intercepted;
eta_r = P_transmitted / P_feed_total;
G_val = D_val * eta_r;
G_dBi = 10*log10(G_val);

fprintf('\n===== Q6 Results =====\n');
fprintf('D_max (uniform):  %.1f dBi\n', D_max_dBi);
fprintf('Directivity D:    %.1f dBi\n', D_dBi);
fprintf('Spillover eff:    %.1f%%\n', eta_spill*100);
fprintf('Reflection eff:   %.1f%%\n', eta_refl*100);
fprintf('Total eta_r:      %.1f%%\n', eta_r*100);
fprintf('Gain G:           %.1f dBi\n', G_dBi);
fprintf('Aperture eff:     %.1f%%\n', D_val/D_max_val*100);
fprintf('======================\n');
