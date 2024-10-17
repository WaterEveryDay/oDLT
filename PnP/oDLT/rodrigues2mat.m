function R = rodrigues2mat(g)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Convert Rodrigues parameters to a rotation matrix.
% Rodrigues parameters are such that g = tan(theta/2) * e,
% where e is the Euler vector (axis of rotation).
%
% References:
% - [1] Markley, F. L., & Crassidis, J. L. (2014). Fundamentals of
%    Spacecraft Attitude Determination and Control (Vol. 33). Springer
%
% Inputs:
% - g (3x1): Rodrigues vector (tangent of half the rotation angle
% times the unit axis of rotation)
%
% Outputs:
% - R (3x3): Rotation matrix


g1 = g(1); g2 = g(2); g3 = g(3);
denom = 1+(g1^2+g2^2+g3^2);

R = 1/denom * ...
    [1 + g1^2 - g2^2 - g3^2, 2*(g1*g2 + g3)        , 2*(g1*g3 - g2);
     2*(g2*g1 - g3)        , 1 - g1^2 + g2^2 - g3^2, 2*(g2*g3 + g1);
     2*(g3*g1 + g2)        , 2*(g3*g2 - g1)        , 1 - g1^2 - g2^2 + g3^2];
end

