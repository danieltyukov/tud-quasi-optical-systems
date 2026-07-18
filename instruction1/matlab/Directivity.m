function [Dir, Prad] = Directivity(E_tot, Theta, dth, dph, er, r)
% Directivity - Compute directivity and total radiated power
%   [Dir, Prad] = Directivity(E_tot, Theta, dth, dph, er, r)

    zeta = 120 * pi / sqrt(er);

    U = (abs(E_tot).^2) / (2 * zeta) * r^2;

    Prad = sum(sum(U .* abs(sin(Theta)) * dth * dph));

    Dir = U / (Prad / (4 * pi));

end
