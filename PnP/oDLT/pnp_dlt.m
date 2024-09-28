function [R, t] = pnp_dlt(X, U, K)
% X = nx3
% U = nx3

% R = rotation from world to camera
% r = position of camera in world

A = build_dlt_system(X, U);

% solve the system by finding smallest eigen vector of A^T A
[~, ~,V] = svd(A, "econ");
h = V(:,end);
H = reshape(h, 3, 4);

% un-calibrate H
H = K \ H;

% transform back to SE3
[R, t] = h_to_se3(H(:));
end
