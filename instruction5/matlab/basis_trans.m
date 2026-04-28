function J = basis_trans(ky, w)
% FT of rect(y/w)/w; sinc(ky*w/2).
    arg = ky * w / 2;
    J   = ones(size(arg));
    nz  = abs(arg) > 1e-12;
    J(nz) = sin(arg(nz)) ./ arg(nz);
end
