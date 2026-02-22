%% Q1: Elementary Electric Source (3 points)

clear; close all; clc;

set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

%% Parameters
c = 3e8;
f = 30e9;
lambda0 = c / f;
k0 = 2*pi / lambda0;
er = 1;

%% Spectral variable
N = 1001;
kx = linspace(0, 3*k0, N);
ky = zeros(size(kx));

%% Compute SGF
sgf = EJ_SGF(er, k0, kx, ky);

%% Q1.1: Plot x-, y-, z-components (real and imaginary parts)

kx_norm = kx / k0;

figure('Color', 'w', 'Position', [100 100 900 700]);

subplot(3,1,1);
plot(kx_norm, real(sgf.Gxx), 'b-', 'LineWidth', 1.5); hold on;
plot(kx_norm, imag(sgf.Gxx), 'r--', 'LineWidth', 1.5);
xline(1, 'k:', 'k_x = k_0', 'LineWidth', 1, 'LabelOrientation', 'horizontal');
xlabel('k_x / k_0');
ylabel('G_{xx}');
title('x-component of SGF (G_{xx})');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

subplot(3,1,2);
plot(kx_norm, real(sgf.Gyx), 'b-', 'LineWidth', 1.5); hold on;
plot(kx_norm, imag(sgf.Gyx), 'r--', 'LineWidth', 1.5);
xline(1, 'k:', 'k_x = k_0', 'LineWidth', 1, 'LabelOrientation', 'horizontal');
xlabel('k_x / k_0');
ylabel('G_{yx}');
title('y-component of SGF (G_{yx})');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

subplot(3,1,3);
plot(kx_norm, real(sgf.Gzx), 'b-', 'LineWidth', 1.5); hold on;
plot(kx_norm, imag(sgf.Gzx), 'r--', 'LineWidth', 1.5);
xline(1, 'k:', 'k_x = k_0', 'LineWidth', 1, 'LabelOrientation', 'horizontal');
xlabel('k_x / k_0');
ylabel('G_{zx}');
title('z-component of SGF (G_{zx})');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

sgtitle('Q1: SGF Components for x-oriented Elementary Source (f = 30 GHz, k_y = 0)');
print('-dpng', '-r150', 'Q1_SGF_components.png');

%% Q1.2: In which region of the spectrum is Gxx imaginary?

fprintf('\n=== Q1.2 Answer ===\n');
fprintf('Gxx is purely imaginary for kx > k0 (evanescent region).\n');
fprintf('In this region kz becomes purely imaginary, making the\n');
fprintf('prefactor -zeta/(2*k*kz) imaginary. Since (k^2-kx^2) is real,\n');
fprintf('Gxx is purely imaginary. For kx < k0 it is purely real.\n');
