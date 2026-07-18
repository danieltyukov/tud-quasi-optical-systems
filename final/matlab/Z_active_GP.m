function Z = Z_active_GP(theta0, phi0, k0, dx, dy, w, l, h, Mmax)
% Active input impedance of an infinite array of x-directed printed dipoles
% over a PEC backing reflector, at scan (theta0, phi0).
%
% Z_act = -1/(dx*dy) * sum_{mx,my} G_xx_eff(kxm, kym) * |I(kxm)|^2 * |Jt(kym)|^2
%
% with G_xx_eff = G_xx_FS * (1 - exp(-j*2*kz*h)) (image theorem above PEC).

    kx0 = k0 * sin(theta0) * cos(phi0);
    ky0 = k0 * sin(theta0) * sin(phi0);

    [mx, my] = meshgrid(-Mmax:Mmax, -Mmax:Mmax);
    kxm = kx0 - 2*pi*mx/dx;
    kym = ky0 - 2*pi*my/dy;

    G   = EJ_SGF_GP(k0, kxm, kym, h);
    Ikx = basis_long(kxm, k0, l);
    Jky = basis_trans(kym, w);

    integrand = G.Gxx .* (Ikx.^2) .* (Jky.^2);
    Z = -sum(integrand(:)) / (dx*dy);
end
