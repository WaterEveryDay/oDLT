function [R, t] = pnp_odlt_vfast_normalized(X, U, K, Sig_uu)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with optimal DLT
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - X (nx3): 3D coordinates of the points in the world frame
% - U (nx3): Pixel coordinates of the corresponding points
% - K (3x3): Camera intrinsic calibration matrix
% - Sig_uu (3x3): Covariance matrix of pixel measurements
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame

if nargin < 4
    Sig_uu = [1, 0, 0; 0, 1, 0; 0, 0, 0];
end

n = size(X, 1);
K_inv = invert_K(K);

U = U*K_inv';
K = eye(3);
K_inv = eye(3);
% normalize the measurements and 3D points for a better behaviour
[X_scaled, T_X, T_X_inv] = normalize_points3D(X);
[U_scaled, T_U, T_U_inv] = normalize_meas(U);

% build the matrix for DLT system
A = build_dlt_system(X_scaled, U_scaled);

% solve for an initial guess using a subset of the system
n_small = min(2*n, min( max(40, floor(n/10)), 100) );
% n_small = 2*n;
idx_small = randperm(2*n, n_small); % floor(linspace(1, 2*n, n_small));

[~,~,V] = svd(A(idx_small, :), 'econ');

h = V(:,end);
H_init = reshape(h, 3, 4);

% compute weights for the DLT system
sig_u = sqrt(Sig_uu(1,1));
scale_3 = H_init(3,:) * [X_scaled, ones(n, 1)]';
q = 1./ (sig_u*scale_3);

% apply weights to the DLT system
A = kron(q', [1; 1]) .* A;

% solve the system by finding smallest eigen vector of A^T A
[~, D, V] = svd(A, 'econ');
h = V(:,end);
inv_cov_h = V * (D.^2) * V';

% bring solution back to regular space
[h, inv_cov_h] = denormalize_kronecker(h, inv_cov_h, K_inv, T_U_inv, T_X);

% solve the Weighted Orthogonal Procrustes Problem
[R, t] = h_to_se3_with_cov_procrustes(h, inv_cov_h);

end