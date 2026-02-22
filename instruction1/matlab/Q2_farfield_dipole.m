%% Q2: Farfield of a Dipole (4 points)

clear; close all; clc;

set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

%% Constants
c = 3e8;
f = 30e9;
lambda0 = c / f;
k0 = 2*pi / lambda0;
er = 1;
zeta = 120 * pi;
R_FF = 1;

%% Q2.1: 1D far-field in main planes
L = lambda0 / 2;
W = lambda0 / 40;

Ntheta = 361;
theta = linspace(0, pi/2, Ntheta);

phi_cuts = [0, pi/4, pi/2];

figure('Color', 'w', 'Position', [100 100 800 500]);
colors = {'b', 'r', 'k'};
styles = {'-', '--', '-.'};

for ip = 1:length(phi_cuts)
    phi = phi_cuts(ip);

    kxs = k0 * sin(theta) * cos(phi);
    kys = k0 * sin(theta) * sin(phi);
    kzs = k0 * cos(theta);

    sgf = EJ_SGF(er, k0, kxs, kys);
    Jx = FTCurrent(k0, er, kxs, kys, L, W);

    [Eth, Eph] = farfield(k0, R_FF, theta, phi*ones(size(theta)), kzs, ...
                          sgf.Gxx, sgf.Gyx, sgf.Gzx, Jx);

    E_tot = sqrt(abs(Eth).^2 + abs(Eph).^2);
    E_tot_dB = 20*log10(E_tot / max(E_tot));

    plot(theta*180/pi, E_tot_dB, [colors{ip} styles{ip}], 'LineWidth', 1.5);
    hold on;
end

xlabel('\theta [deg]');
ylabel('Normalized |E^{far}| [dB]');
title('Q2.1: Far Field in Main Planes (L = \lambda_0/2, W = \lambda_0/40, f = 30 GHz)');
legend('\phi = 0°', '\phi = 45°', '\phi = 90°', 'Location', 'southwest');
grid on;
ylim([-30 0]);
xlim([0 90]);
print('-dpng', '-r150', 'Q2_1_farfield_1D.png');

%% Q2.2: UV representation
Nth = 181;
Nph = 361;
theta_2d = linspace(0, pi/2, Nth);
phi_2d = linspace(0, 2*pi, Nph);
[TH2, PH2] = meshgrid(theta_2d, phi_2d);

kxs2 = k0 * sin(TH2) .* cos(PH2);
kys2 = k0 * sin(TH2) .* sin(PH2);
kzs2 = k0 * cos(TH2);

sgf2 = EJ_SGF(er, k0, kxs2, kys2);
Jx2 = FTCurrent(k0, er, kxs2, kys2, L, W);

[Eth2, Eph2] = farfield(k0, R_FF, TH2, PH2, kzs2, ...
                         sgf2.Gxx, sgf2.Gyx, sgf2.Gzx, Jx2);

E_tot2 = sqrt(abs(Eth2).^2 + abs(Eph2).^2);
E_tot2_dB = 20*log10(E_tot2 / max(E_tot2(:)));

U = sin(TH2) .* cos(PH2);
V = sin(TH2) .* sin(PH2);

figure('Color', 'w', 'Position', [100 100 700 600]);
surf(U, V, E_tot2_dB, 'EdgeColor', 'none');
view(2);
colormap('jet');
caxis([-10 0]);
colorbar;
xlabel('U = sin\theta cos\phi');
ylabel('V = sin\theta sin\phi');
title('Q2.2: Far Field in UV Plane (10 dB dynamic range)');
axis equal;
xlim([-1 1]);
ylim([-1 1]);
print('-dpng', '-r150', 'Q2_2_farfield_UV.png');

%% Q2.3: Broadside directivity vs frequency
L_fixed = 5e-3;
W_fixed = 0.25e-3;

freq = linspace(5e9, 60e9, 200);
Dir_broadside_dB = zeros(size(freq));

Nth_dir = 91;
Nph_dir = 181;
theta_dir = linspace(0, pi/2, Nth_dir);
phi_dir = linspace(0, 2*pi, Nph_dir);
dth_dir = theta_dir(2) - theta_dir(1);
dph_dir = phi_dir(2) - phi_dir(1);
[TH_dir, PH_dir] = meshgrid(theta_dir, phi_dir);

for ifreq = 1:length(freq)
    f_i = freq(ifreq);
    k_i = 2*pi*f_i / c;

    kxs_i = k_i * sin(TH_dir) .* cos(PH_dir);
    kys_i = k_i * sin(TH_dir) .* sin(PH_dir);
    kzs_i = k_i * cos(TH_dir);

    sgf_i = EJ_SGF(er, k_i, kxs_i, kys_i);
    Jx_i = FTCurrent(k_i, er, kxs_i, kys_i, L_fixed, W_fixed);

    [Eth_i, Eph_i] = farfield(k_i, R_FF, TH_dir, PH_dir, kzs_i, ...
                               sgf_i.Gxx, sgf_i.Gyx, sgf_i.Gzx, Jx_i);

    E_tot_i = sqrt(abs(Eth_i).^2 + abs(Eph_i).^2);

    U_i = (abs(E_tot_i).^2) / (2 * zeta) * R_FF^2;

    Prad_upper = sum(sum(U_i .* sin(TH_dir) * dth_dir * dph_dir));
    Prad_total = 2 * Prad_upper;

    U_broadside = U_i(1, 1);
    Dir_bs = U_broadside / (Prad_total / (4*pi));

    Dir_broadside_dB(ifreq) = 10*log10(max(Dir_bs, 1e-20));
end

figure('Color', 'w', 'Position', [100 100 800 500]);
plot(freq/1e9, Dir_broadside_dB, 'b-', 'LineWidth', 1.5);
xlabel('Frequency [GHz]');
ylabel('Directivity at broadside [dBi]');
title('Q2.3: Broadside Directivity vs Frequency (L = 5 mm, W = 0.25 mm)');
grid on;
hold on;
yline(2.15, 'r--', 'Half-wave dipole (2.15 dBi)', 'LineWidth', 1);
f_res = c / (2 * L_fixed);
xline(f_res/1e9, 'k:', ['f_{res} = ' num2str(f_res/1e9) ' GHz'], ...
       'LineWidth', 1, 'LabelOrientation', 'horizontal');
print('-dpng', '-r150', 'Q2_3_directivity_vs_freq.png');
