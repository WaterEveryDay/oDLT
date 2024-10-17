function [R, t] = pnp_dlt_normalized(X, U, K)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with the normalized DLT
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - X (nx3): 3D coordinates of the points in the world frame
% - U (nx3): Pixel coordinates of the corresponding points
% - K (3x3): Camera intrinsic calibration matrix
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame

n = size(X, 1);
K_inv = invert_K(K);

X = [X, ones(n, 1)];

% center and scale points
[U_scaled, ~, T_U_inv] = normalize_meas(U);
[X_scaled, T_X] = normalize_points3D(X);

% build the system
A = build_dlt_system(X_scaled, U_scaled);

% solve the system by finding smallest eigen vector of A^T A
[~,~,V] = svd(A, "econ");
h = V(:,end);
H = reshape(h, 3, 4);

% de-normalize and un-calibrate the camera matrix
H = K_inv * (T_U_inv * H * T_X);

% transform back to SE3
[R, t] = h_to_se3(H(:));

end
