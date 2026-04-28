%% Q2: Grating-lobe diagrams for the two lattices.
% Plots Floquet-mode circles (radius k0, centers at (2*pi*mx/dx, 2*pi*my/dy))
% in normalized (kx0/k0, ky0/k0) space and overlays the visible region (unit
% circle) plus the scan trajectories for E, H and D plane.
clear; close all; clc;
set_plot_defaults();

%% Parameters
freq   = 10e9;
c      = 3e8;
lambda = c/freq;
k0     = 2*pi/lambda;

cases = struct( ...
    'name', {'15 mm', '20 mm'}, ...
    'dx',   {15e-3, 20e-3}, ...
    'dy',   {15e-3, 20e-3});

planes = struct( ...
    'name',  {'E-plane (\phi = 0\circ)', 'H-plane (\phi = 90\circ)', 'D-plane (\phi = 45\circ)'}, ...
    'phi',   {0, pi/2, pi/4}, ...
    'color', {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]}, ...
    'style', {'-', '--', '-.'});

%% Plot
fig = figure('Color','w','Position',[80 80 1500 780]);
set(fig, 'InvertHardcopy', 'off');

th_circ = linspace(0, 2*pi, 361);
cx = cos(th_circ); cy = sin(th_circ);

xlims = [-2.6 2.6];
ylims = [-2.6 2.6];

for ci = 1:numel(cases)
    dx = cases(ci).dx;
    dy = cases(ci).dy;
    rx = lambda/dx;
    ry = lambda/dy;

    subplot(1, numel(cases), ci); hold on; axis equal; grid on; box on;

    % --- Higher-order Floquet circles ---
    Mplot = 2;
    for mx = -Mplot:Mplot
        for my = -Mplot:Mplot
            if mx == 0 && my == 0, continue; end
            cxm = mx * rx;
            cym = my * ry;
            % Only draw the part that is inside the plot window
            inside = (cxm + cx >= xlims(1)) & (cxm + cx <= xlims(2)) & ...
                     (cym + cy >= ylims(1)) & (cym + cy <= ylims(2));
            if any(inside)
                plot(cxm + cx, cym + cy, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
                if cxm >= xlims(1) && cxm <= xlims(2) && cym >= ylims(1) && cym <= ylims(2)
                    plot(cxm, cym, 'k.', 'MarkerSize', 7);
                    text(cxm + 0.07, cym + 0.10, sprintf('(%d,%d)', mx, my), ...
                         'FontSize', 8, 'Color', [0.3 0.3 0.3]);
                end
            end
        end
    end

    % --- Visible region (fundamental mode), drawn on top ---
    fill(cx, cy, [0.85 0.92 1.0], 'EdgeColor', [0 0.3 0.7], 'LineWidth', 2, ...
         'FaceAlpha', 0.55);
    plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
    text(0.07, -0.12, '(0,0)', 'FontSize', 9, 'FontWeight', 'bold');

    % --- Scan trajectories: solid arrow segments inside visible region ---
    th_scan = linspace(-1, 1, 201);
    h_lines = gobjects(1, numel(planes));
    for pi_ = 1:numel(planes)
        u = th_scan * cos(planes(pi_).phi);
        v = th_scan * sin(planes(pi_).phi);
        h_lines(pi_) = plot(u, v, planes(pi_).style, ...
             'LineWidth', 2.4, 'Color', planes(pi_).color);
    end

    % --- Mark GL-entry points if any ---
    sinE = rx - 1;  sinH = ry - 1;
    if sinE >= 0 && sinE <= 1
        plot([+sinE, -sinE], [0, 0], 'rp', 'MarkerSize', 16, ...
             'MarkerFaceColor', [1 0.3 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        if abs(sinE - 1) < 0.01
            text(+sinE - 0.4, -0.25, sprintf('GL@\\theta=±%.0f°', asind(sinE)), 'FontSize', 9);
        else
            text(+sinE + 0.05, -0.28, sprintf('GL@\\theta=±%.0f° (E)', asind(sinE)), 'FontSize', 9);
        end
    end
    if sinH >= 0 && sinH <= 1
        plot([0, 0], [+sinH, -sinH], 'rp', 'MarkerSize', 16, ...
             'MarkerFaceColor', [1 0.3 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        if abs(sinH - 1) < 0.01
            % already labelled above (axes coincide visually)
        else
            text(0.07, +sinH + 0.18, sprintf('GL@\\theta=±%.0f° (H)', asind(sinH)), 'FontSize', 9);
        end
    end

    xlabel('k_{x0}/k_0  =  sin\theta_0 cos\phi_0');
    ylabel('k_{y0}/k_0  =  sin\theta_0 sin\phi_0');
    title(sprintf('d_x = d_y = %s   (\\lambda/d = %.3f)', cases(ci).name, rx));
    xlim(xlims); ylim(ylims);
    xticks(-2:1:2); yticks(-2:1:2);

    if ci == 1
        legend(h_lines, {planes.name}, 'Location', 'southoutside', 'NumColumns', 3);
    end
end

sgtitle('Grating-lobe diagrams (10 GHz): Floquet-mode circles and scan trajectories', ...
        'FontSize', 12, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile('..','figures','Q2_grating_lobe_diagram.png'));

%% Print the threshold angles for the report
fprintf('\nGrating-lobe entry analysis at f = 10 GHz, lambda = %g mm:\n', lambda*1e3);
for ci = 1:numel(cases)
    rx = lambda/cases(ci).dx; ry = lambda/cases(ci).dy;
    fprintf('  %s lattice (lambda/d = %.4f):\n', cases(ci).name, rx);
    % E-plane GL entry: sin(th) = lambda/dx - 1 (requires <=1, i.e. dx >= lambda/2)
    sinE = rx - 1;
    if sinE > 1
        fprintf('    E-plane GL entry: never (lambda/d - 1 = %.3f > 1)\n', sinE);
    else
        fprintf('    E-plane GL entry: theta = %6.2f deg  (mode (-1,0) or (+1,0))\n', asind(sinE));
    end
    % H-plane GL entry: sin(th) = lambda/dy - 1
    sinH = ry - 1;
    if sinH > 1
        fprintf('    H-plane GL entry: never (lambda/d - 1 = %.3f > 1)\n', sinH);
    else
        fprintf('    H-plane GL entry: theta = %6.2f deg  (mode (0,-1) or (0,+1))\n', asind(sinH));
    end
    % D-plane: nearest mode (-1,0) gives equation
    %   (u/sqrt2 + rx)^2 + (u/sqrt2)^2 = 1   with u = sin(theta)
    %   u^2 + sqrt(2)*rx*u + rx^2 - 1 = 0
    a_ = 1; b_ = sqrt(2)*rx; cD = rx^2 - 1;
    disc = b_^2 - 4*a_*cD;
    if disc < 0
        fprintf('    D-plane GL entry: never (mode (-1,0) miss: disc=%.3f<0)\n', disc);
    else
        roots_ = roots([a_ b_ cD]);
        uD = roots_(abs(roots_) <= 1);
        if isempty(uD)
            fprintf('    D-plane GL entry: never (no root in [-1,1])\n');
        else
            fprintf('    D-plane GL entry: theta = %6.2f deg  (mode (-1,0))\n', asind(max(uD)));
        end
    end
end
