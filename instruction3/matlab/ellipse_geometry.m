function lens = ellipse_geometry(er, D, lambda0)
% Elliptical lens geometry from permittivity and diameter.

    lens.e = 1/sqrt(er);
    lens.b = D/2;
    lens.a = lens.b / sqrt(1 - lens.e^2);
    lens.c = lens.e * lens.a;
    lens.h = lens.a * (1 + lens.e);
    lens.theta_max = atan2(lens.b, lens.e*lens.a);
    lens.D = D;
    lens.lambda0 = lambda0;
    lens.er = er;
end
