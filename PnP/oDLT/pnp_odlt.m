function [R, t] = ... %hbar_uncalibrated, cov_hbar_uncalibrated] = ...
    pnp_opt_dlt(X, U, K, Sig_uu)
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
[X_scaled, T_X, T_X_inv] = normalize_points3D(X);
[U_scaled, T_U, T_U_inv] = normalize_meas(U);

% build the matrix for DLT system
A = build_dlt_system(X_scaled, U_scaled);

% solve for an initial guess using a subset of the system
n_small = min(2*n, min( max(30, floor(n/10)), 100) );

% rng(1)
[~,~,V] = svd(A(randperm(2*n, n_small), :), "econ");

h = V(:,end);
H_init = reshape(h, 3, 4);

% de-normalize and un-calibrate the camera matrix
H = K_inv * (T_U_inv * H_init * T_X);

% transform back to SE3
[R, t] = h_to_se3(H(:));
H = K * [R, t];
H_init = T_U * H * T_X_inv;

% compute weights for the DLT system
sig_u = T_U(1,1)*sqrt(Sig_uu(1,1));
X_proj = H_init * [X_scaled, ones(n, 1)]';
scale_3 = X_proj(3,:);
inv_scale_norm_sig = 1./ (scale_3 * sig_u);

% apply weights to the DLT system
A = kron(inv_scale_norm_sig', [1; 1]) .* A;

% solve the system by finding smallest eigen vector of A^T A
[~, Sig, V] = svd(A, 'econ');
h = V(:,end);
% cov_h = inv(A'*A);
cov_h = V * diag(1 ./ diag(Sig.^2)) * V';

% compute vec(inv(K) * inv(T_U) * H * T_X)
% and its covariance
[h, cov_h] = denormalize_kronecker(h, cov_h, K_inv, T_U_inv, T_X);

% solve the Weighted Orthogonal Procrustes Problem
[R, t] = h_to_se3_with_cov_procrustes(h, cov_h);

end
    


