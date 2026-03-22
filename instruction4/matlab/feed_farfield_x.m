function [E_theta, E_phi] = feed_farfield_x(k0, a_feed, theta, phi)
% feed_farfield_x - Far-field pattern of circular aperture feed (x-polarized)
%   [E_theta, E_phi] = feed_farfield_x(k0, a_feed, theta, phi)
%   The feed is a uniform circular aperture with x-directed electric current.
%   k0     : free-space wavenumber
%   a_feed : feed aperture radius
%   theta  : elevation angles (can be array)
%   phi    : azimuth angles (same size as theta)
%
%   The far field is computed using the spectral Green's function approach:
%   E = j*kz * G^ej * J * exp(-jkR)/(2*pi*R)
%   For x-directed current, use the x-column of G^ej: Gxx, Gyx, Gzx

    R_FF = 1;  % unit distance for pattern (cancels in normalization)
    er = 1;    % free space

    kx = k0 * sin(theta) .* cos(phi);
    ky = k0 * sin(theta) .* sin(phi);
    kz = k0 * cos(theta);
    k_perp = k0 * sin(theta);

    % Spectral Green's function
    sgf = EJ_SGF(er, k0, kx, ky);

    % FT of uniform x-directed circular current
    Jx = FTCircular(a_feed, k_perp);

    % Far field using x-column of SGF: Gxx, Gyx, Gzx
    [E_theta, E_phi] = farfield_ej(k0, R_FF, theta, phi, kz, ...
                                    sgf.Gxx, sgf.Gyx, sgf.Gzx, Jx);

end
