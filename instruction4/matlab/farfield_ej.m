function [Eth, Eph] = farfield_ej(k0, R_FF, TH, PH, kz, G1, G2, G3, J)
% Far-field E_theta, E_phi via stationary phase.
% G1,G2,G3: any column of SGF (x-current: Gxx,Gyx,Gzx; y-current: Gxy,Gyy,Gzy)

    phase_factor = exp(-1j * k0 * R_FF) ./ (2 * pi * R_FF);

    jkz_G1 = 1j .* kz .* G1;
    jkz_G2 = 1j .* kz .* G2;
    jkz_G3 = 1j .* kz .* G3;

    bad = isnan(jkz_G1) | isnan(jkz_G2) | isnan(jkz_G3) | ...
          isinf(jkz_G1) | isinf(jkz_G2) | isinf(jkz_G3);
    jkz_G1(bad) = 0;
    jkz_G2(bad) = 0;
    jkz_G3(bad) = 0;

    Ex = jkz_G1 .* J .* phase_factor;
    Ey = jkz_G2 .* J .* phase_factor;
    Ez = jkz_G3 .* J .* phase_factor;

    Eth = Ex .* cos(TH) .* cos(PH) + Ey .* cos(TH) .* sin(PH) - Ez .* sin(TH);
    Eph = -Ex .* sin(PH) + Ey .* cos(PH);
end
