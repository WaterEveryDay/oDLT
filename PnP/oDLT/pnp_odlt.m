function [R, t] = ... %hbar_uncalibrated, cov_hbar_uncalibrated] = ...
    pnp_odlt(X, U, K, Sig_uu)
% X = 3xn: 3D points in world
% U = 3xn: 2D pixel measurements
% K = 3x3: Calibration matrix
% Sig_uu = 3x3: Covariance of pixel measurements

% R = 3 x 3 rotation from world to camera
% r = 3 x 1 position of camera in world

% arguments
%     X(:,3) double
%     U(:,3) double
%     K(3,3) double
%     Sig_uu(3,3) double = [1, 0, 0; 0, 1, 0; 0, 0, 0];
% end
if nargin < 4
    Sig_uu = [1, 0, 0; 0, 1, 0; 0, 0, 0];
end

n = size(X, 1);
K_inv = invert_K(K);

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
q = 1./ (scale_3 * sig_u);

% apply weights to the DLT system
A = kron(q', [1; 1]) .* A;

% solve the system by finding smallest eigen vector of A^T A
[~, Sig, V] = svd(A, 'econ');
h = V(:,end);
inv_cov_h = V * (Sig.^2) * V';

% compute vec(inv(K) * inv(T_U) * H * T_X)
% and its covariance
[h, inv_cov_h] = denormalize_kronecker(h, inv_cov_h, K_inv, T_U_inv, T_X);

% solve the Weighted Orthogonal Procrustes Problem
[R, t] = h_to_se3_with_cov_procrustes(h, inv_cov_h);
end