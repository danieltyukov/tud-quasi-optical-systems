function Z = Z_FSS(theta0, phi0, k0, dx, dy, w, l, Mmax)
% FSS Galerkin self-impedance, single sinusoidal basis function, free space.
% Same Floquet sum as Z_active (instruction 5): the source type does not
% change Z, only the right-hand side v.
    kx0 = k0 * sin(theta0) * cos(phi0);
    ky0 = k0 * sin(theta0) * sin(phi0);

    [mx, my] = meshgrid(-Mmax:Mmax, -Mmax:Mmax);
    kxm = kx0 - 2*pi*mx/dx;
    kym = ky0 - 2*pi*my/dy;

    G   = EJ_SGF(1, k0, kxm, kym);
    Ikx = basis_long(kxm, k0, l);
    Jky = basis_trans(kym, w);

    integrand = G.Gxx .* (Ikx.^2) .* (Jky.^2);
    Z = -sum(integrand(:)) / (dx*dy);
end
