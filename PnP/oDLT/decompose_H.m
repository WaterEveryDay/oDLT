function [R, t, K] = decompose_H(H)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% This function extracts the rotation matrix, translation vector, 
% and intrinsic calibration matrix from an unscaled camera matrix.
%
%
% Inputs:
% - H (3x4): Unscaled camera matrix
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame
% - K (3x3): Camera intrinsic calibration matrix


H_inv = inv(H(1:3,1:3));

Rot = [-1, 0, 0; 0, -1, 0; 0, 0, 1];
% recover the inverse rotation and inverse of K
[R_transpose, K_inv] = qr(H_inv);

% recover the calibration matrix
K = inv(K_inv);

% scale such that the last element of the calibration matrix is 1
scale = K(3,3);
K = K*Rot/scale;

% recover extrinsics
R = Rot*R_transpose';
r = - H_inv * H(:,4);
t = - R * r;
end