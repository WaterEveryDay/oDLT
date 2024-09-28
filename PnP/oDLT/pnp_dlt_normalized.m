function [R, t] = pnp_dlt_normalized(X, U, K)
% X = n x 3
% U = n x 3

% R = rotation from world to camera
% r = position of camera in world

n = size(X, 1);
K_inv = invert_K(K);

X = [X, ones(n, 1)];
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
