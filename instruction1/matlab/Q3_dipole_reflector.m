%% Q3: Dipole with a Backing Reflector (3 points)

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

%% Dipole parameters
L = lambda0 / 2;
W = lambda0 / 40;

%% Q3.1: Far field with PEC reflector at h = 7.5 mm

h = 7.5e-3;

Ntheta = 361;
theta = linspace(0, pi/2, Ntheta);

phi_cuts = [0, pi/4, pi/2];
phi_labels = {'0', '45', '90'};

figure('Color', 'w', 'Position', [100 100 900 600]);

for ip = 1:length(phi_cuts)
    phi = phi_cuts(ip);

    kxs = k0 * sin(theta) * cos(phi);
    kys = k0 * sin(theta) * sin(phi);
    kzs = k0 * cos(theta);

    sgf = EJ_SGF(er, k0, kxs, kys);
    Jx = FTCurrent(k0, er, kxs, kys, L, W);

    % Array factor from image theorem
    AF = 2j * sin(kzs * h);

    [Eth_ref, Eph_ref] = farfield(k0, R_FF, theta, phi*ones(size(theta)), kzs, ...
                                   sgf.Gxx, sgf.Gyx, sgf.Gzx, Jx .* AF);
    E_tot_ref = sqrt(abs(Eth_ref).^2 + abs(Eph_ref).^2);

    [Eth_fs, Eph_fs] = farfield(k0, R_FF, theta, phi*ones(size(theta)), kzs, ...
                                 sgf.Gxx, sgf.Gyx, sgf.Gzx, Jx);
    E_tot_fs = sqrt(abs(Eth_fs).^2 + abs(Eph_fs).^2);

    ref_max = max(E_tot_fs);
    E_ref_dB = 20*log10(E_tot_ref / ref_max);
    E_fs_dB = 20*log10(E_tot_fs / ref_max);

    subplot(1,3,ip);
    plot(theta*180/pi, E_fs_dB, 'b-', 'LineWidth', 1.5); hold on;
    plot(theta*180/pi, E_ref_dB, 'r--', 'LineWidth', 1.5);
    xlabel('\theta [deg]');
    ylabel('|E^{far}| [dB]');
    title(['\phi = ' phi_labels{ip} '°']);
    legend('Free space', 'With PEC reflector', 'Location', 'southwest');
    grid on;
    ylim([-30 10]);
    xlim([0 90]);
end

sgtitle('Q3.1: Farfield Comparison - Free Space vs PEC Reflector (h = 7.5 mm)');
print('-dpng', '-r150', 'Q3_1_farfield_comparison.png');

%% Q3.3: Directivity at broadside vs distance h

Nh = 200;
h_range = linspace(lambda0/10, 2*lambda0, Nh);

Nth_dir = 181;
Nph_dir = 361;
theta_dir = linspace(0, pi/2, Nth_dir);
phi_dir = linspace(0, 2*pi, Nph_dir);
dth_dir = theta_dir(2) - theta_dir(1);
dph_dir = phi_dir(2) - phi_dir(1);
[TH_dir, PH_dir] = meshgrid(theta_dir, phi_dir);

kxs_dir = k0 * sin(TH_dir) .* cos(PH_dir);
kys_dir = k0 * sin(TH_dir) .* sin(PH_dir);
kzs_dir = k0 * cos(TH_dir);

sgf_dir = EJ_SGF(er, k0, kxs_dir, kys_dir);
Jx_dir = FTCurrent(k0, er, kxs_dir, kys_dir, L, W);

Dir_broadside = zeros(size(h_range));

for ih = 1:length(h_range)
    h_i = h_range(ih);

    AF_i = 2j * sin(kzs_dir * h_i);

    [Eth_i, Eph_i] = farfield(k0, R_FF, TH_dir, PH_dir, kzs_dir, ...
                               sgf_dir.Gxx, sgf_dir.Gyx, sgf_dir.Gzx, Jx_dir .* AF_i);

    E_tot_i = sqrt(abs(Eth_i).^2 + abs(Eph_i).^2);

    [Dir_i, ~] = Directivity(E_tot_i, TH_dir, dth_dir, dph_dir, er, R_FF);
    Dir_broadside(ih) = Dir_i(1, 1);
end

figure('Color', 'w', 'Position', [100 100 800 500]);
plot(h_range/lambda0, Dir_broadside, 'b-', 'LineWidth', 1.5);
xlabel('h / \lambda_0');
ylabel('Directivity at broadside (linear)');
title('Q3.3: Broadside Directivity vs Distance from Ground Plane');
grid on;

hold on;
for n = 1:4
    h_null = n * lambda0 / 2;
    if h_null <= 2*lambda0
        xline(h_null/lambda0, 'r:', ['h = ' num2str(n) '\lambda_0/2'], ...
               'LineWidth', 1, 'LabelOrientation', 'horizontal');
    end
end

print('-dpng', '-r150', 'Q3_3_directivity_vs_h.png');

%% Q3.3: Verify Prad = 0 for h = 0
AF_h0 = 2j * sin(kzs_dir * 0);
[Eth_h0, Eph_h0] = farfield(k0, R_FF, TH_dir, PH_dir, kzs_dir, ...
                              sgf_dir.Gxx, sgf_dir.Gyx, sgf_dir.Gzx, Jx_dir .* AF_h0);
E_tot_h0 = sqrt(abs(Eth_h0).^2 + abs(Eph_h0).^2);
Prad_h0 = sum(sum(abs(E_tot_h0).^2 / (2*zeta) * R_FF^2 .* abs(sin(TH_dir)) * dth_dir * dph_dir));
fprintf('Prad(h=0) = %.6e W\n', Prad_h0);
