function v = v_FSS(theta0, phi0, V_TM, V_TE, k0, l, w)
% Galerkin "voltage" for plane-wave illumination. Only the fundamental
% Floquet mode contributes (the incident field is a single plane wave).
    kx0 = k0 * sin(theta0) * cos(phi0);
    ky0 = k0 * sin(theta0) * sin(phi0);

    Ekx_inc = cos(theta0)*cos(phi0)*V_TM - sin(phi0)*V_TE;
    B = basis_long(kx0, k0, l) .* basis_trans(ky0, w);

    v = Ekx_inc .* B;
end
