function Voc = compute_Voc(Va_th, Va_ph, VGO_th, VGO_ph, theta, dth, dphi, zeta)
% compute_Voc - Open circuit voltage of antenna in focal plane (reception analysis)
%
%   Voc = compute_Voc(Va_th, Va_ph, VGO_th, VGO_ph, theta, dth, dphi, zeta)
%
%   Va_th, Va_ph   : theta/phi components of antenna Tx far field (angular, 2D grid)
%                    These are the "voltage" patterns V_a^tx = E_a^tx * r * exp(jkr)
%   VGO_th, VGO_ph : theta/phi components of GO field on FO sphere (angular, 2D grid)
%                    These are V_r^GO = e_r^GO * R * exp(jkR)
%   theta          : 2D meshgrid of theta values
%   dth            : theta step
%   dphi           : phi step
%   zeta           : free-space impedance (120*pi)
%
%   The open-circuit voltage is (from lecture notes):
%     Voc = (2/zeta) * integral{ V_a^tx . V_r^GO * sin(theta) dtheta dphi }
%
%   The dot product is computed component-wise in theta,phi basis.

    integrand = (Va_th .* VGO_th + Va_ph .* VGO_ph) .* sin(theta);

    Voc = (2 / zeta) * sum(integrand(:)) * dth * dphi;

end
