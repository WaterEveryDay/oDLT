function [R, t, K] = pnp_odlt_uncalibrated(X, U, Sig_uu)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with optimal DLT in the case of an uncalibrated camera
%
% References
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - X (nx3): 3D coordinates of the points in the world frame
% - U (nx3): Pixel coordinates of the corresponding points
% - Sig_uu (3x3): Covariance matrix of pixel measurements
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame
% - K (3x3): Camera intrinsic calibration matrix

if nargin < 4
    Sig_uu = [1, 0, 0; 0, 1, 0; 0, 0, 0];
end

n = size(X, 1);

% normalize the measurements and 3D points for a better behaviour
[X_scaled, T_X, ~] = normalize_points3D(X);
[U_scaled, T_U, T_U_inv] = normalize_meas(U);

% build the matrix for DLT system
A = build_dlt_system(X_scaled, U_scaled);

% solve for an initial guess using a subset of the system
n_small = min(2*n, min( max(30, floor(n/10)), 100) );
% n_small = 2*n;
[~,~,V] = svd(A(randperm(2*n, n_small), :), "econ");

h = V(:,end);
H_init = reshape(h, 3, 4);

% compute weights for the DLT system
sig_u = T_U(1,1)*sqrt(Sig_uu(1,1));
X_proj = H_init * [X_scaled, ones(n, 1)]';
scale_3 = X_proj(3,:);
inv_scale_norm_sig = 1./ (scale_3 * sig_u);

% apply weights to the DLT system
A = kron(inv_scale_norm_sig', [1; 1]) .* A;

% solve the system by finding smallest eigen vector of A^T A
[~, ~, V] = svd(A, 'econ');
h = V(:,end);

% de-normalize and un-calibrate the camera matrix
H = reshape(h, 3, 4);
H = T_U_inv * H * T_X;
[R, t, K] = decompose_H(H);
end
    


