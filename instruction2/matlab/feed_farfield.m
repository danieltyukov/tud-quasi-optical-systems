function [E_theta, E_phi] = feed_farfield(k0, a_feed, theta, phi)
% feed_farfield - Far-field pattern of circular aperture feed (y-polarized)
%   [E_theta, E_phi] = feed_farfield(k0, a_feed, theta, phi)
%   k0     : free-space wavenumber
%   a_feed : feed aperture radius
%   theta  : elevation angles (can be array)
%   phi    : azimuth angles (same size as theta)

    R_FF = 1;
    er = 1;

    kx = k0 * sin(theta) .* cos(phi);
    ky = k0 * sin(theta) .* sin(phi);
    kz = k0 * cos(theta);
    k_perp = k0 * sin(theta);

    % Spectral Green's function
    sgf = EJ_SGF(er, k0, kx, ky);

    % FT of uniform y-directed circular current
    Jy = FTCircular(a_feed, k_perp);

    % Far field using y-column of SGF: Gxy, Gyy, Gzy
    [E_theta, E_phi] = farfield_ej(k0, R_FF, theta, phi, kz, ...
                                    sgf.Gxy, sgf.Gyy, sgf.Gzy, Jy);

end
