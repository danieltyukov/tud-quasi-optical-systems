function theta_i = lens_incidence_angle(e, theta)
% Local incidence angle on elliptical lens surface.

    cos_i = (1 - e*cos(theta)) ./ sqrt(1 + e.^2 - 2*e*cos(theta));
    cos_i = min(max(cos_i, -1), 1);
    theta_i = acos(cos_i);
end
