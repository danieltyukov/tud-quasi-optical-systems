function Jx = FTCurrent(k0, er, kx, ky, l, w)
% FTCurrent - Fourier Transform of x-directed dipole current distribution
%   Jx = FTCurrent(k0, er, kx, ky, l, w)

    keq = k0 * sqrt(er);

    % T(ky) = sinc(ky*w/2)
    ky_arg = ky * w / 2;
    T = ones(size(ky_arg));
    idx = abs(ky_arg) > 1e-15;
    T(idx) = sin(ky_arg(idx)) ./ ky_arg(idx);

    % L(kx)
    numerator = 2 * keq * (cos(kx * l / 2) - cos(keq * l / 2));
    denominator = (keq^2 - kx.^2) * sin(keq * l / 2);

    L_kx = zeros(size(kx));
    singular = abs(keq^2 - kx.^2) < 1e-10;
    L_kx(~singular) = numerator(~singular) ./ denominator(~singular);

    if any(singular)
        kx_near = keq * (1 + 1e-8);
        num_near = 2 * keq * (cos(kx_near * l / 2) - cos(keq * l / 2));
        den_near = (keq^2 - kx_near^2) * sin(keq * l / 2);
        L_kx(singular) = num_near / den_near;
    end

    Jx = L_kx .* T;

end
