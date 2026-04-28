function I = basis_long(kx, k0, l)
% FT of sinusoidal current sin(k0*(l/2-|x|))/sin(k0*l/2) on a dipole of length l.
    s   = sin(k0*l/2);
    num = 2*k0 .* (cos(kx*l/2) - cos(k0*l/2));
    den = (k0^2 - kx.^2) .* s;
    I = num ./ den;

    % Removable 0/0 at |kx|=k0; L'Hopital limit is l/2.
    near_pole = abs(k0^2 - kx.^2) < 1e-10 * k0^2;
    I(near_pole) = l/2;
end
