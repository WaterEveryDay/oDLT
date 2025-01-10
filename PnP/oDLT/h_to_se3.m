function [R, t] = h_to_se3(h)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% This function extracts the rotation matrix and translation vector 
% from an unscaled and unorthogonalized camera matrix.
%
% Inputs:
% - h (12x1): Unscaled and unorthogonalized camera matrix 
%   H = [R, -R*r], where R is the rotation matrix and r is the 
%   camera position vector in the world frame.
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector


% extract the rotation matrix and compute the scaling
R_init = reshape(h(1:9), 3, 3);
detR = det(R_init);
scale = sign(detR) * abs(detR)^(1/3);

% scale h
h = h/scale;
R_init = reshape(h(1:9), 3, 3);
[U, ~, V] = svd(R_init);

% retrieve SO3 matrix and translation
R = U *diag([1, 1, det(U*V')])* V';
t = h(10:12);
end