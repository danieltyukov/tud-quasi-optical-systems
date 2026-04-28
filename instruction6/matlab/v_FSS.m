function v = v_FSS(theta0, phi0, V_TM, V_TE, k0, l, w)
% v_FSS - Galerkin "voltage" for an FSS of dipole strips under plane-wave
% illumination.  For an incident plane wave with TM amplitude V_TM and TE
% amplitude V_TE, only the fundamental Floquet mode contributes:
%
%   v = ( cos(theta)*cos(phi)*V_TM - sin(phi)*V_TE ) * B(-kx0, -ky0)
%
% where B(kx,ky) = I(kx)*Jt(ky) is the FT of the basis function.  I and Jt
% are even in their argument, so B(-kx0,-ky0) = I(kx0)*Jt(ky0).
%
% Reference: notes section 3.4 "Field Definitions and Incident Field
% Decomposition".  Sign convention: the incident E-field projected on x is
% E_inc_x = cos(theta)*cos(phi)*V_TM - sin(phi)*V_TE.

    kx0 = k0 * sin(theta0) * cos(phi0);
    ky0 = k0 * sin(theta0) * sin(phi0);

    Ekx_inc = cos(theta0)*cos(phi0)*V_TM - sin(phi0)*V_TE;
    B = basis_long(kx0, k0, l) .* basis_trans(ky0, w);

    v = Ekx_inc .* B;
end
