function F = FTCircular(a, k_perp)
% FTCircular - 2D Fourier transform of a uniform circular aperture
%   F = FTCircular(a, k_perp)
%   a      : radius of circular aperture
%   k_perp : transverse wavenumber = sqrt(kx^2 + ky^2)
%   Output : pi*a^2 * 2*J1(k_perp*a) / (k_perp*a), with limit pi*a^2 at k_perp=0

    F = zeros(size(k_perp));
    idx = abs(k_perp) > 1e-15;
    arg = k_perp(idx) * a;
    F(idx) = pi * a^2 * 2 * besselj(1, arg) ./ arg;
    F(~idx) = pi * a^2;

end
