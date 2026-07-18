%% Q3+Q4: Windowed FSS pattern, 10x10 cells, normal vs oblique TM incidence.
clear; close all; clc;
set_plot_defaults();

freq    = 10e9;
c       = 3e8;
lambda  = c/freq;
k0      = 2*pi/lambda;
dx      = 20e-3;
dy      = 20e-3;
w       = 1e-3;
l       = 15e-3;
Mmax    = 30;
Nx      = 10;
Ny      = 10;
phi0    = 0;
V_TM    = 1; V_TE = 0;

theta_in = asind(lambda/dx - 1);

theta_obs_deg = linspace(-89.9, 89.9, 3601);
theta_obs     = deg2rad(theta_obs_deg);

incidence = struct( ...
    'name',  {'\theta_0 = 0° (normal)', '\theta_0 = \theta_{in} = 30° (GL onset)'}, ...
    'tag',   {'normal', 'oblique'}, ...
    'theta0', {0, deg2rad(theta_in)});

patterns = cell(numel(incidence),1);
peak_dB  = zeros(numel(incidence),1);
GL_dB    = nan(numel(incidence),1);
GL_th    = nan(numel(incidence),1);
iBF_inc  = zeros(numel(incidence),1);
Z_inc    = zeros(numel(incidence),1);
GL_info  = struct();

for ii = 1:numel(incidence)
    th0  = incidence(ii).theta0;
    kx0  = k0*sin(th0)*cos(phi0);
    ky0  = k0*sin(th0)*sin(phi0);

    Z_inc(ii)   = Z_FSS(th0, phi0, k0, dx, dy, w, l, Mmax);
    v_in        = v_FSS(th0, phi0, V_TM, V_TE, k0, l, w);
    iBF_inc(ii) = v_in / Z_inc(ii);

    kx = k0*sin(theta_obs);
    ky = zeros(size(kx));
    kz = k0*cos(theta_obs);

    G   = EJ_SGF(1, k0, kx, ky);
    Ikx = basis_long(kx, k0, l);
    Jky = basis_trans(ky, w);

    % Closed-form Dirichlet kernel; explicit limit at integer multiples of pi.
    psi_x = (kx - kx0) * dx / 2;
    AFx   = ones(size(psi_x)) * Nx;
    nz    = abs(sin(psi_x)) > 1e-12;
    AFx(nz) = sin(Nx*psi_x(nz)) ./ sin(psi_x(nz)) .* exp(1j*(Nx-1)*psi_x(nz));

    psi_y = (ky - ky0) * dy / 2;
    AFy   = ones(size(psi_y)) * Ny;
    nz    = abs(sin(psi_y)) > 1e-12;
    AFy(nz) = sin(Ny*psi_y(nz)) ./ sin(psi_y(nz)) .* exp(1j*(Ny-1)*psi_y(nz));

    E = (1j*kz) .* G.Gxx .* Ikx .* Jky .* AFx .* AFy * iBF_inc(ii);

    P = abs(E).^2;
    P_dB = 10*log10(P);
    P_dB = P_dB - max(P_dB);
    patterns{ii} = struct('th_deg', theta_obs_deg, 'P_dB', P_dB, 'E', E);

    [Pmax, imax] = max(abs(E));
    peak_dB(ii)  = 20*log10(Pmax);
    fprintf('Incidence %s:  Z = %.2f + j(%.2f) Ohm,  i_BF = %.3e,  peak at theta = %+.2f deg\n', ...
            incidence(ii).tag, real(Z_inc(ii)), imag(Z_inc(ii)), abs(iBF_inc(ii)), ...
            theta_obs_deg(imax));

    if ii == 2
        s_gl  = sin(th0) - lambda/dx;
        th_gl_theory = NaN;
        if abs(s_gl) <= 1
            th_gl_theory = asind(s_gl);
        end

        if ~isnan(th_gl_theory)
            [~, iGLth] = min(abs(theta_obs_deg - th_gl_theory));
            P_at_GLtheory = P(iGLth);
        else
            iGLth = NaN; P_at_GLtheory = NaN;
        end

        mask_far = theta_obs_deg < -45;
        [PglFar, iFar] = max(P .* mask_far);
        GL_dB(ii) = 10*log10(PglFar / Pmax^2);
        GL_th(ii) = theta_obs_deg(iFar);

        mask_near    = abs(theta_obs_deg - theta_obs_deg(imax)) > 4 & ...
                       abs(theta_obs_deg - theta_obs_deg(imax)) < 12;
        [Psl, iSL]   = max(P .* mask_near);
        SL1_dB       = 10*log10(Psl / Pmax^2);
        SL1_th       = theta_obs_deg(iSL);

        fprintf('  Theoretical GL: sin(theta_GL) = %.3f -> theta_GL = %.2f deg\n', ...
                s_gl, th_gl_theory);
        fprintf('  Pattern level near GL (theta = %+.2f deg): %.2f dB below main lobe\n', ...
                theta_obs_deg(iGLth), 10*log10(P_at_GLtheory/Pmax^2));
        fprintf('  Highest peak in GL region (theta < -45 deg): theta = %+.2f deg, level = %.2f dB\n', ...
                GL_th(ii), GL_dB(ii));
        fprintf('  First sidelobe near main beam: theta = %+.2f deg, level = %.2f dB\n', ...
                SL1_th, SL1_dB);

        GL_info = struct('th_gl_theory', th_gl_theory, ...
                         'P_at_GLtheory_dB', 10*log10(P_at_GLtheory/Pmax^2), ...
                         'GL_far_th', GL_th(ii), 'GL_far_dB', GL_dB(ii), ...
                         'SL1_th', SL1_th, 'SL1_dB', SL1_dB);
    end
end

%% Patterns, two incidence cases
fig = figure('Color','w','Position',[80 80 1200 480]);
set(fig, 'InvertHardcopy', 'off');

for ii = 1:numel(incidence)
    subplot(1,2,ii); hold on; grid on;
    plot(patterns{ii}.th_deg, patterns{ii}.P_dB, 'LineWidth', 1.6, ...
         'Color', [0 0.45 0.74]);
    xlabel('\theta_{obs}  [deg]'); ylabel('|E(\theta)|^2 / |E_{peak}|^2  [dB]');
    title(incidence(ii).name);
    xlim([-90 90]); ylim([-80 1]); xticks(-90:30:90);

    th0_deg = rad2deg(incidence(ii).theta0);
    xline(th0_deg, 'k:', 'LineWidth', 1, 'Alpha', 0.6);

    if ii == 2
        s_gl = sin(deg2rad(theta_in)) - lambda/dx;
        if abs(s_gl) <= 1
            th_gl = asind(s_gl);
            xline(th_gl, 'r:', 'LineWidth', 1.2);
            text(th_gl+1.0, -8, sprintf('GL at \\theta = %+.0f°', th_gl), ...
                 'FontSize', 9, 'Color', [0.7 0 0]);
        end
        if isfield(GL_info,'th_gl_theory') && ~isnan(GL_info.th_gl_theory)
            plot(GL_info.th_gl_theory, GL_info.P_at_GLtheory_dB, 'rp', ...
                 'MarkerSize', 14, 'MarkerFaceColor', [1 0.3 0], ...
                 'MarkerEdgeColor', 'k');
            text(GL_info.th_gl_theory+2, GL_info.P_at_GLtheory_dB+4, ...
                 sprintf('%.1f dB', GL_info.P_at_GLtheory_dB), ...
                 'FontSize', 9, 'Color', [0.7 0 0]);
        end
    end
end

sgtitle({'Windowed FSS scattered-pattern, 10×10 unit cells, TM, \phi=0 plane', ...
         sprintf('d_x = d_y = %g mm,   l = %g mm,   f = %g GHz   (\\lambda/d = %.2f)', ...
                 dx*1e3, l*1e3, freq/1e9, lambda/dx)}, ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q3_windowed_patterns.png'));

%% Overlay: dB + linear common-normalisation
fig2 = figure('Color','w','Position',[80 80 1100 420]);
set(fig2, 'InvertHardcopy', 'off');
subplot(1,2,1); hold on; grid on;
plot(patterns{1}.th_deg, patterns{1}.P_dB, '-',  'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(patterns{2}.th_deg, patterns{2}.P_dB, '--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
xlabel('\theta_{obs}  [deg]'); ylabel('|E|^2 / |E_{peak}|^2  [dB]');
title('Overlay (each normalised to its own peak)');
legend({incidence.name}, 'Location', 'south');
xlim([-90 90]); ylim([-60 1]); xticks(-90:30:90);

subplot(1,2,2); hold on; grid on;
P0 = abs(patterns{1}.E).^2;
P1 = abs(patterns{2}.E).^2;
gmax = max([max(P0), max(P1)]);
plot(patterns{1}.th_deg, P0/gmax, '-',  'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
plot(patterns{2}.th_deg, P1/gmax, '--', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
xlabel('\theta_{obs}  [deg]'); ylabel('|E|^2  (normalised, linear)');
title('Linear scale (common normalisation)');
xlim([-90 90]); ylim([0 1.05]); xticks(-90:30:90);

sgtitle('Pattern comparison: normal vs oblique incidence (TM, \phi = 0)', ...
        'FontSize', 11, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q3_windowed_overlay.png'));

save('Q3_pattern_results.mat', 'theta_obs_deg','patterns','peak_dB','GL_dB','GL_th', ...
     'iBF_inc','Z_inc','Nx','Ny','dx','dy','w','l','freq','lambda','theta_in');

fprintf('\nFigures saved.\n');
