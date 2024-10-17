function [h, cov_h_inv] = denormalize_kronecker(h, cov_h_inv, K_inv, T_U_inv, T_X)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Denormalize the camera matrix and its inverse covariance, along with
% de-clamping the camera calibration.
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - h (12x1): vectorized camera matrix in the similarity space
% - cov_h_inv (12x12): inverse covariance of h in the similarity space
% - K_inv (3x3): inverse camera calibration matrix
% - T_U_inv (3x3): inverse similarity transformation on the pixel measurements
% - T_X (4x4): similarity transformation on the 3D points
%
% Outputs:
% - h (12x1):  vectorized camera matrix in the regular space
% - cov_h_inv (12x12): inverse covariance of h in the regular space

% Step 1: Create the transformation matrix
M = kron(T_X', K_inv * T_U_inv);

% Step 2: Apply the transformation to h
h = M * h;

M_inv = inv(M);
% Step 3: Transform the covariance matrix
cov_h_inv = M_inv' * cov_h_inv * M_inv;
end

