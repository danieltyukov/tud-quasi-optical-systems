function [Eth, Eph] = farfield(k0, R_FF, TH, PH, kz, Gxx, Gyx, Gzx, Jx)
% farfield - Far-field E_theta, E_phi via stationary phase evaluation
%   [Eth, Eph] = farfield(k0, R_FF, TH, PH, kz, Gxx, Gyx, Gzx, Jx)

    phase_factor = exp(-1j * k0 * R_FF) ./ (2 * pi * R_FF);

    % j*kz * G: multiply kz into the passed SGF components.
    % At theta=pi/2 kz=0 and Gxx/Gyx diverge (1/kz prefactor), but their
    % product with kz is finite.  Replace any resulting NaN with the
    % analytical limit: j*kz*G_ij = -j*zeta/(2k) * matrix_element.
    jkz_Gxx = 1j .* kz .* Gxx;
    jkz_Gyx = 1j .* kz .* Gyx;
    jkz_Gzx = 1j .* kz .* Gzx;

    bad = isnan(jkz_Gxx) | isnan(jkz_Gyx) | isnan(jkz_Gzx);
    if any(bad(:))
        zeta = 120 * pi;
        kxs = k0 * sin(TH(bad)) .* cos(PH(bad));
        kys = k0 * sin(TH(bad)) .* sin(PH(bad));
        common = -1j * zeta / (2 * k0);
        jkz_Gxx(bad) = common .* (k0^2 - kxs.^2);
        jkz_Gyx(bad) = common .* (-kxs .* kys);
        jkz_Gzx(bad) = common .* (-kxs .* kz(bad));
    end

    Ex = jkz_Gxx .* Jx .* phase_factor;
    Ey = jkz_Gyx .* Jx .* phase_factor;
    Ez = jkz_Gzx .* Jx .* phase_factor;

    Eth = Ex .* cos(TH) .* cos(PH) + Ey .* cos(TH) .* sin(PH) - Ez .* sin(TH);
    Eph = -Ex .* sin(PH) + Ey .* cos(PH);

end
