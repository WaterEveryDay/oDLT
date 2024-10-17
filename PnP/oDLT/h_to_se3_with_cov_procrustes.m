function [R, t] = h_to_se3_with_cov_procrustes(h, cov_h_inv)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% This function extracts the rotation matrix and translation vector 
% from an unscaled and unorthogonalized camera matrix.
% Accounts for covariance information
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - h (12x1): Unscaled and unorthogonalized camera matrix 
%   H = [R, -R*r], where R is the rotation matrix and r is the 
%   translation vector.
% - cov_h_inv (12x12): the inverse covariance of h
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame

% extract the rotation matrix and compute the scaling
R_init = reshape(h(1:9), 3, 3);
detR = det(R_init);
scale = sign(detR) * abs(detR)^(1/3);

% scale h
h = h/scale;
R_init = reshape(h(1:9), 3, 3);

if rank(cov_h_inv) < 12
    [U, ~, V] = svd(R_init);
    R = U * diag([1, 1, det(U * V')]) * V';
else
    W = reshape(diag(cov_h_inv(1:9, 1:9)), 3, 3);
    R = weighted_orthogonal_procrustes_gn(R_init, W);
end

t = h(10:12);
end