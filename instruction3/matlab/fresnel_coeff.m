function [tau_perp, tau_par, Gamma_perp, Gamma_par] = fresnel_coeff(er, theta_i)
% Fresnel coefficients for dielectric-to-air interface.

    zeta_0 = 120*pi;
    zeta_d = zeta_0/sqrt(er);

    sin_t = sqrt(er)*sin(theta_i);
    cos_t = sqrt(1 - sin_t.^2);
    cos_i = cos(theta_i);

    Gamma_perp = (zeta_0*cos_i - zeta_d*cos_t) ./ (zeta_0*cos_i + zeta_d*cos_t);
    tau_perp   = (2*zeta_0*cos_i) ./ (zeta_0*cos_i + zeta_d*cos_t);

    Gamma_par = (zeta_0*cos_t - zeta_d*cos_i) ./ (zeta_0*cos_t + zeta_d*cos_i);
    tau_par   = (2*zeta_0*cos_i) ./ (zeta_0*cos_t + zeta_d*cos_i);
end
