function R = rotVecToMat(phi)
%   Converts an Euler rotation vector to an Euler axis/angle
%   pair (angle in radians).

theta = norm(phi);

if theta < 1e-12
    e = [0;0;0];
else
    e = phi./theta;
end

R = cos(theta)*eye(3) - sin(theta)*cross_mat(e) + (1 - cos(theta))*(e*e');
end

