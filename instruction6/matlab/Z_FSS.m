function Z = Z_FSS(theta0, phi0, k0, dx, dy, w, l, Mmax)
% Z_FSS - Self-impedance of an infinite FSS of dipole strips, single basis
% function (Galerkin), free space.
%   Z = -1/(dx*dy) * sum_{mx,my} G_xx^{ej}(kxm, kym) * |B(kxm, kym)|^2
%   with kxm = kx0 - 2*pi*mx/dx, kym = ky0 - 2*pi*my/dy and the basis
%   B(kx,ky) = I(kx)*Jt(ky) (real), so |B|^2 = I(kx)^2 * Jt(ky)^2.
%
% This is mathematically identical to the active-array input impedance
% (instruction 5, Z_active.m); the FSS uses the same Floquet sum because
% the periodicity and Galerkin projection do not depend on whether the
% excitation is a feed voltage or a plane wave.

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
