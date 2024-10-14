function [R, t] = ... %hbar_uncalibrated, cov_hbar_uncalibrated] = ...
    pnp_odlt_full(X, U, K, Sig_uu, Sig_pp_all)
% X = 3xn: 3D points in world
% U = 3xn: 2D pixel measurements
% K = 3x3: Calibration matrix
% Sig_uu = 3x3: Covariance of pixel measurements
% Sig_pp_all = 3x3xn: Covariance of 3D points

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

if nargin < 5
end

n = size(X, 1);
K_inv = invert_K(K);

% normalize the measurements and 3D points for a better behaviour
[X_scaled, T_X, ~] = normalize_points3D(X);
[U_scaled, T_U, T_U_inv] = normalize_meas(U);

% build the matrix for DLT system
A = build_dlt_system_full(X_scaled, U_scaled);

% solve for an initial guess using a subset of the system
n_small = min(2*n, min( max(30, floor(n/10)), 100) );
[~,~,V] = svd(A(randperm(2*n, n_small), :), "econ");

h = V(:,end);
H_init = reshape(h, 3, 4);

% compute weights for the DLT system
Sig_uu_scaled = T_U * Sig_uu * T_U';
scale_3 = H_init(end,:) * [X_scaled, ones(n, 1)]';

cross_u_all = zeros(3,3,n);
cross_u_all(1,2,:) = -U_scaled(:,3);
cross_u_all(2,1,:) = U_scaled(:,3);

cross_u_all(1,3,:) = U_scaled(:,2);
cross_u_all(3,1,:) = -U_scaled(:,2);

cross_u_all(2,3,:) = -U_scaled(:,1);
cross_u_all(3,2,:) = U_scaled(:,1);

% compute the combination of covariance of 3D point and pixel measurement
% on the residual
Sig_combined = reshape(scale_3.^2, 1, 1, n) .* Sig_uu_scaled;
if nargin >= 5
    Sig_xx_scaled = T_X .* Sig_pp_all .* T_X';
    Sig_combined = Sig_combined + H_init(1:3,1:3)*Sig_pp_scaled_all*H_init(1:3,1:3)';
end
Sig_eps_all = - pagemtimes(pagemtimes(cross_u_all, Sig_combined), cross_u_all);
Sig_eps_pinv_all = pagepinv(Sig_eps_all);

A_all = pagetranspose(reshape(A', 12, 3, n));
ASigA = pagemtimes(pagemtimes(pagetranspose(A_all), Sig_eps_pinv_all), A_all);
sumA = sum(ASigA, 3);

% solve the system by finding smallest eigen vector of A^T A
[~, Sig, V] = svd(sumA, 'econ');
h = V(:,end);
inv_cov_h = V * (Sig) * V';

% compute vec(inv(K) * inv(T_U) * H * T_X)
% and its covariance
[h, inv_cov_h] = denormalize_kronecker(h, inv_cov_h, K_inv, T_U_inv, T_X);

% solve the Weighted Orthogonal Procrustes Problem
[R, t] = h_to_se3_with_cov_procrustes(h, inv_cov_h);
end