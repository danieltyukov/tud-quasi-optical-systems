function [ej_SGF] = EJ_SGF(er, k, kx, ky)
% EJ_SGF - Spectral dyadic Green's function for electric currents
%   [ej_SGF] = EJ_SGF(er, k, kx, ky)

    zeta = 120 * pi / sqrt(er);

    kz = -1j * sqrt(-(k^2 - kx.^2 - ky.^2));

    prefactor = -zeta ./ (2 * k * kz);

    ej_SGF.Gxx = prefactor .* (k^2 - kx.^2);
    ej_SGF.Gxy = prefactor .* (-kx .* ky);
    ej_SGF.Gxz = prefactor .* (-kx .* kz);
    ej_SGF.Gyx = prefactor .* (-kx .* ky);
    ej_SGF.Gyy = prefactor .* (k^2 - ky.^2);
    ej_SGF.Gyz = prefactor .* (-ky .* kz);
    ej_SGF.Gzx = prefactor .* (-kx .* kz);
    ej_SGF.Gzy = prefactor .* (-ky .* kz);
    ej_SGF.Gzz = prefactor .* (k^2 - kz.^2);

end
