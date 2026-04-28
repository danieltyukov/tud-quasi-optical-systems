function [G] = EJ_SGF_GP(k0, kx, ky, h)
% Spectral Green's function for an x-directed electric current at height h
% above an infinite PEC ground plane, evaluated at the source plane (z = z' = h).
%
% Image theorem: PEC at z=0 induces an opposite-sign image of the horizontal
% electric current at z = -h. The total field at z = h is therefore
%   G_eff = G_FS(z'=h) - G_FS(z'=-h)
% In the spectral domain, that becomes
%   G_eff = G_FS * (1 - exp(-j*2*kz*h))
% with the kz branch chosen so evanescent modes decay (Im{kz} < 0).
    zeta = 120 * pi;
    kz   = -1j * sqrt(-(k0^2 - kx.^2 - ky.^2));   % propagating: real > 0; evanescent: -j*|.|
    prefactor = -zeta ./ (2 * k0 * kz);

    Gxx_FS = prefactor .* (k0^2 - kx.^2);

    % Reflector factor: 0 at h=0 (image cancels source), 2 at h=lambda/4 broadside.
    refl = 1 - exp(-1j * 2 * kz * h);

    G.Gxx = Gxx_FS .* refl;
    G.kz  = kz;
    G.refl_factor = refl;
end
