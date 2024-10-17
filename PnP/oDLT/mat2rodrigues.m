function g = mat2rodrigues(R)
% Sebastien Henry
%
% Convert rotation matrix to Rodrigues parameters.
% Rodrigues parameters are such that g = tan(theta/2) * e,
% where e is the Euler vector (axis of rotation).
%
% References:
% - [1] Markley, F. L., & Crassidis, J. L. (2014). Fundamentals of
%    Spacecraft Attitude Determination and Control (Vol. 33). Springer
%
% Inputs:
% - R (3x3): Rotation matrix
%
% Outputs:
% - g (3x1): Rodrigues vector (tangent of half the rotation angle
% times the unit axis of rotation)


theta = real(acos((trace(R)-1)/2));
eps = 1e-12;

if abs(theta) < eps
    % if angle is very small, then rodrigues is zero
    g = [0;0;0];
    return
elseif abs(theta - pi) < 2*eps
    % if angle is close to pi, then singularity
    % perturb with minor angle to avoid singularity
    theta = theta + 2*eps;
end

e = 1/(2*sin(theta))*[R(2,3)-R(3,2), R(3,1)-R(1,3), R(1,2)-R(2,1)]';
g = tan(theta/2) * e;

end

