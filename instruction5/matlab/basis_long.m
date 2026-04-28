function I = basis_long(kx, k0, l)
% basis_long - FT of the longitudinal sinusoidal current on a dipole of length l.
%   i(x) = sin(k0*(l/2 - |x|))/sin(k0*l/2)  for |x| <= l/2
%   I(kx) = 2*k0*(cos(kx*l/2) - cos(k0*l/2)) / ((k0^2 - kx^2)*sin(k0*l/2))
%   At kx -> ±k0 the expression has a removable 0/0; the L'Hôpital limit is l/2.

    s   = sin(k0*l/2);
    num = 2*k0 .* (cos(kx*l/2) - cos(k0*l/2));
    den = (k0^2 - kx.^2) .* s;

    I = num ./ den;

    % Removable pole at |kx| = k0 -> limit = l/2
    near_pole = abs(k0^2 - kx.^2) < 1e-10 * k0^2;
    I(near_pole) = l/2;
end
