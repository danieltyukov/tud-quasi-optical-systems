function [E_theta, E_phi] = feed_farfield_x(k0, a_feed, theta, phi)
% Far-field pattern of circular aperture feed (x-polarized uniform current)

    R_FF = 1;
    er = 1;

    kx = k0 * sin(theta) .* cos(phi);
    ky = k0 * sin(theta) .* sin(phi);
    kz = k0 * cos(theta);
    k_perp = k0 * sin(theta);

    sgf = EJ_SGF(er, k0, kx, ky);

    % FT of uniform x-directed circular current
    Jx = FTCircular(a_feed, k_perp);

    % x-column of SGF: Gxx, Gyx, Gzx
    [E_theta, E_phi] = farfield_ej(k0, R_FF, theta, phi, kz, ...
                                    sgf.Gxx, sgf.Gyx, sgf.Gzx, Jx);

end
