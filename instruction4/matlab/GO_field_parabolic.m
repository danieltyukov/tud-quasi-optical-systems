function [e_GO_theta, e_GO_phi] = GO_field_parabolic(E0_PW, theta, phi, theta_0)
% GO_field_parabolic - GO electric field on the FO sphere for a parabolic reflector
%   illuminated by an x-polarized broadside plane wave.
%
%   [e_GO_theta, e_GO_phi] = GO_field_parabolic(E0_PW, theta, phi, theta_0)
%
%   E0_PW   : amplitude of incident plane wave [V/m]
%   theta   : elevation angles on FO sphere (0 to theta_0)
%   phi     : azimuth angles on FO sphere (0 to 2*pi)
%   theta_0 : reflector half-angle (rim angle)
%
%   The GO field for an x-polarized broadside PW on a parabolic reflector:
%     e_r^GO = -S_par * (cos(phi)*theta_hat - sin(phi)*phi_hat) * E0_PW
%   where S_par = 2/(1+cos(theta)) is the spreading factor.
%
%   Output: theta and phi components of the GO field (spectral/angular domain).
%   The field is set to zero outside theta_0.

    % Spreading factor for parabolic reflector
    S_par = 2 ./ (1 + cos(theta));

    % GO field components (from lecture notes eq. for broadside x-pol PW)
    e_GO_theta = -S_par .* cos(phi) .* E0_PW;
    e_GO_phi   =  S_par .* sin(phi) .* E0_PW;

    % Zero field outside the reflector rim
    mask = theta <= theta_0;
    e_GO_theta = e_GO_theta .* mask;
    e_GO_phi   = e_GO_phi .* mask;

end
