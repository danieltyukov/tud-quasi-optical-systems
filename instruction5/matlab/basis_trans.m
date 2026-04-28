function J = basis_trans(ky, w)
% basis_trans - FT of the transverse rect current of width w (normalized to width).
%   j_t(y) = rect_{[-w/2, w/2]}(y) / w  ->  J_t(ky) = sinc(ky*w/2)  (unnormalized)

    arg = ky * w / 2;
    J   = ones(size(arg));
    nz  = abs(arg) > 1e-12;
    J(nz) = sin(arg(nz)) ./ arg(nz);
end
