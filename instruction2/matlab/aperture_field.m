function [E_a_x, E_a_y] = aperture_field(k0, f_focal, a_feed, theta_p, phi_p)
% GO aperture field on parabolic reflector.
% Projects feed far field onto aperture plane with amplitude taper.

    [E_th_feed, E_ph_feed] = feed_farfield(k0, a_feed, theta_p, phi_p);

    % Amplitude taper: (1+cos(theta'))/(2f) from 1/r spreading
    amplitude_taper = (1 + cos(theta_p)) ./ (2 * f_focal);

    % Project spherical (theta,phi) components onto Cartesian (x,y)
    E_a_x = amplitude_taper .* (E_th_feed .* cos(theta_p) .* cos(phi_p) - E_ph_feed .* sin(phi_p));
    E_a_y = amplitude_taper .* (E_th_feed .* cos(theta_p) .* sin(phi_p) + E_ph_feed .* cos(phi_p));
end
