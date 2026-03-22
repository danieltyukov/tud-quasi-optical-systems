function [e_GO_theta, e_GO_phi] = GO_field_off_broadside(E0_PW, k0, F, theta, phi, theta_0, theta_i, phi_i)
% GO_field_off_broadside - GO field on FO sphere for off-broadside plane wave
%
%   [e_GO_theta, e_GO_phi] = GO_field_off_broadside(E0_PW, k0, F, theta, phi, theta_0, theta_i, phi_i)
%
%   E0_PW   : amplitude of incident x-polarized plane wave [V/m]
%   k0      : free-space wavenumber
%   F       : focal length of reflector (= FO sphere radius R)
%   theta   : elevation angles on FO sphere
%   phi     : azimuth angles on FO sphere
%   theta_0 : reflector rim angle
%   theta_i : incidence angle of the plane wave (from broadside)
%   phi_i   : azimuth angle of incidence
%
%   For slightly off-broadside incidence, the GO field is the broadside
%   field multiplied by a phase due to the transverse path difference:
%
%     e_r^GO(off) = e_r^GO(bs) * exp(-j * Delta_k_rho . rho_refl)
%
%   where rho_refl = 2F*tan(theta/2) is the radial position on the
%   reflector aperture, and Delta_k_rho is the transverse wavenumber shift.
%
%   This can be decomposed as:
%     Linear phase:  exp(-j * Delta_k_rho . rho_FO)  [beam steering]
%     Coma phase:    exp(-j * Delta_k_rho . (rho_refl - rho_FO) )  [aberration]
%   where rho_FO = F*sin(theta) = R*sin(theta).

    % Handle broadside case quickly
    if theta_i < 1e-15
        [e_GO_theta, e_GO_phi] = GO_field_parabolic(E0_PW, theta, phi, theta_0);
        return;
    end

    % Broadside GO field
    S_par = 2 ./ (1 + cos(theta));
    e_GO_theta_bs = -S_par .* cos(phi) .* E0_PW;
    e_GO_phi_bs   =  S_par .* sin(phi) .* E0_PW;

    % Transverse wave vector of the incident PW
    dk_x = k0 * sin(theta_i) * cos(phi_i);
    dk_y = k0 * sin(theta_i) * sin(phi_i);

    % Position on reflector aperture: rho' = 2F*tan(theta/2)
    % In Cartesian: x' = rho'*cos(phi), y' = rho'*sin(phi)
    rho_refl = 2 * F * tan(theta/2);
    x_refl = rho_refl .* cos(phi);
    y_refl = rho_refl .* sin(phi);

    % Total phase shift from off-broadside incidence
    % The incident PW has a phase exp(-j*(dk_x*x + dk_y*y)) on the aperture
    total_phase = exp(-1j * (dk_x * x_refl + dk_y * y_refl));

    % Apply phase to broadside GO field
    e_GO_theta = e_GO_theta_bs .* total_phase;
    e_GO_phi   = e_GO_phi_bs   .* total_phase;

    % Zero field outside the reflector rim
    mask = theta <= theta_0;
    e_GO_theta = e_GO_theta .* mask;
    e_GO_phi   = e_GO_phi .* mask;

end
