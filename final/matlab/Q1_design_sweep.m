%% Q1 design sweep: choose (h, l, dx, Z_L) to meet |Gamma_act| < -10 dB
%% over Ka-band 27.5-31 GHz at broadside.
%%
%% We fix w = 0.4 mm (thin printed strip) and dx = dy = 4.8 mm
%% (just below lambda_min/2 at f_max = 31 GHz, where lambda_min/2 = 4.84 mm)
%% to be grating-lobe-free up to large scan angles in the entire band.
%%
%% Strategy: sweep (h, l) on a coarse grid, compute the worst-case |Gamma|
%% across the band for each Z_L candidate, and pick the design that gives
%% the largest in-band margin below -10 dB.

clear; close all; clc;
set_plot_defaults();

c0    = 3e8;
f_lo  = 27.5e9; f_hi = 31e9; f_c = (f_lo+f_hi)/2;
lam_c = c0/f_c;
lam_lo = c0/f_lo;
lam_hi = c0/f_hi;

w     = 0.4e-3;
dx    = 4.8e-3;  dy = 4.8e-3;
Mmax  = 25;

% Frequency grid in band (used for the band-edge diagnostic)
freqs_band = linspace(f_lo, f_hi, 31);

% --- Coarse sweep over (h, l)
h_vals = linspace(0.10, 0.30, 21) * lam_c;     % around lambda/4
l_vals = linspace(0.30, 0.48, 19) * lam_c;     % half-wave at center is 0.5*lam_c, cell size limits this

% For each (h,l) we sweep Z_L over a range and pick the best
ZL_vals = 50:10:300;

worst_GammaDB = nan(numel(h_vals), numel(l_vals), numel(ZL_vals));

fprintf('Sweeping (h, l, Z_L) for broadside band performance...\n');
for ih = 1:numel(h_vals)
    for il = 1:numel(l_vals)
        h_now = h_vals(ih);
        l_now = l_vals(il);
        if l_now > 0.95*dx, continue; end       % keep dipole inside the cell
        Z_band = zeros(size(freqs_band));
        for kf = 1:numel(freqs_band)
            k0 = 2*pi*freqs_band(kf)/c0;
            Z_band(kf) = Z_active_GP(0, 0, k0, dx, dy, w, l_now, h_now, Mmax);
        end
        for iz = 1:numel(ZL_vals)
            G = (Z_band - ZL_vals(iz)) ./ (Z_band + ZL_vals(iz));
            worst_GammaDB(ih, il, iz) = 20*log10(max(abs(G)));
        end
    end
    fprintf('  h = %5.3f mm done (%d/%d)\n', h_vals(ih)*1e3, ih, numel(h_vals));
end

% --- Pick the design with smallest worst-case |Gamma|
[wmin, idx] = min(worst_GammaDB(:));
[ih_best, il_best, iz_best] = ind2sub(size(worst_GammaDB), idx);
h_best  = h_vals(ih_best);
l_best  = l_vals(il_best);
ZL_best = ZL_vals(iz_best);

fprintf('\nBest broadside design:\n');
fprintf('  h   = %.3f mm  (= %.3f * lambda_c)\n', h_best*1e3, h_best/lam_c);
fprintf('  l   = %.3f mm  (= %.3f * lambda_c)\n', l_best*1e3, l_best/lam_c);
fprintf('  w   = %.3f mm,   dx = dy = %.3f mm\n', w*1e3, dx*1e3);
fprintf('  Z_L = %.0f Ohm\n', ZL_best);
fprintf('  worst |Gamma| in band = %.2f dB\n', wmin);

% --- Recompute on a dense band grid for plotting
freqs = linspace(0.85*f_lo, 1.10*f_hi, 200);
Z_f   = zeros(size(freqs));
for kf = 1:numel(freqs)
    k0 = 2*pi*freqs(kf)/c0;
    Z_f(kf) = Z_active_GP(0, 0, k0, dx, dy, w, l_best, h_best, Mmax);
end
Gamma_f = (Z_f - ZL_best) ./ (Z_f + ZL_best);

% Compare to NO reflector (set h -> infinity) and to a simple lambda/4 default
Z_f_FS = zeros(size(freqs));
for kf = 1:numel(freqs)
    k0 = 2*pi*freqs(kf)/c0;
    [mx, my] = meshgrid(-Mmax:Mmax, -Mmax:Mmax);
    G_fs = EJ_SGF(1, k0, -2*pi*mx/dx, -2*pi*my/dy);
    Ikx  = basis_long(-2*pi*mx/dx, k0, l_best);
    Jky  = basis_trans(-2*pi*my/dy, w);
    Z_f_FS(kf) = -sum(sum(G_fs.Gxx .* (Ikx.^2) .* (Jky.^2)))/(dx*dy);
end
Gamma_f_FS = (Z_f_FS - ZL_best)./(Z_f_FS + ZL_best);

%% Plot 1 : |Gamma| vs frequency at broadside
fig = figure('Color','w','Position',[80 80 1100 700]);
set(fig, 'InvertHardcopy', 'off');

subplot(2,2,1); hold on; grid on;
plot(freqs/1e9, real(Z_f), 'b-', 'LineWidth', 1.7);
plot(freqs/1e9, imag(Z_f), 'r--', 'LineWidth', 1.7);
yline(ZL_best, 'k:', sprintf('Z_L = %.0f \\Omega', ZL_best), 'LabelHorizontalAlignment', 'left');
xline(f_lo/1e9, 'k-.','LineWidth', 0.7); xline(f_hi/1e9, 'k-.','LineWidth', 0.7);
xlabel('frequency [GHz]'); ylabel('Z_{in,act} [\Omega]');
title('Active input impedance at broadside');
legend({'Re\{Z\}','Im\{Z\}'}, 'Location', 'best');
xlim([min(freqs) max(freqs)]/1e9);

subplot(2,2,2); hold on; grid on;
plot(freqs/1e9, 20*log10(abs(Gamma_f)),    'b-',  'LineWidth', 1.8);
plot(freqs/1e9, 20*log10(abs(Gamma_f_FS)), 'r--', 'LineWidth', 1.4);
yline(-10, 'k:', '-10 dB');
xline(f_lo/1e9, 'k-.','LineWidth', 0.7); xline(f_hi/1e9, 'k-.','LineWidth', 0.7);
xlabel('frequency [GHz]'); ylabel('|\Gamma_{act}| [dB]');
title('Active reflection coefficient at broadside');
legend({'with backing reflector','free space (no reflector)'}, 'Location', 'best');
ylim([-30 0]); xlim([min(freqs) max(freqs)]/1e9);

% Map of (h,l) sweep at the chosen Z_L
subplot(2,2,3);
[Lg, Hg] = meshgrid(l_vals*1e3, h_vals*1e3);
contourf(Lg, Hg, worst_GammaDB(:,:,iz_best), -25:1:0, 'LineColor', 'none');
colorbar; colormap(parula);
hold on;
plot(l_best*1e3, h_best*1e3, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
[Cf, hC] = contour(Lg, Hg, worst_GammaDB(:,:,iz_best), [-10 -10], 'k-', 'LineWidth', 1.5);
clabel(Cf, hC, 'FontSize', 9, 'Color', 'k');
xlabel('dipole length l [mm]'); ylabel('reflector height h [mm]');
title(sprintf('worst |\\Gamma| in band [dB], Z_L = %d \\Omega', ZL_best));

% Sweep over Z_L at the best (h, l)
subplot(2,2,4); hold on; grid on;
worst_at_best = squeeze(worst_GammaDB(ih_best, il_best, :));
plot(ZL_vals, worst_at_best, 'b-o', 'LineWidth', 1.6, 'MarkerSize', 4);
plot(ZL_best, wmin, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
yline(-10, 'k:', '-10 dB');
xlabel('Z_L [\Omega]'); ylabel('worst |\Gamma_{act}| in band [dB]');
title('Z_L optimisation at best (h, l)');

sgtitle({sprintf('Q1: Ka-band Satcom dipole array — design at broadside (f = %.2f-%.2f GHz)', f_lo/1e9, f_hi/1e9), ...
    sprintf('h = %.2f mm,  l = %.2f mm,  w = %.2f mm,  d_x = d_y = %.2f mm,  Z_L = %d \\Omega', ...
            h_best*1e3, l_best*1e3, w*1e3, dx*1e3, ZL_best)}, ...
    'FontSize', 11, 'FontWeight', 'bold');
print('-dpng','-r150', fullfile('..','figures','Q1_design_sweep.png'));

save('Q1_design.mat', ...
     'h_best','l_best','w','dx','dy','ZL_best','Mmax', ...
     'freqs','Z_f','Gamma_f','f_lo','f_hi','f_c','lam_c', ...
     'worst_GammaDB','h_vals','l_vals','ZL_vals','wmin');

fprintf('\nFigure saved to ../figures/Q1_design_sweep.png\n');
