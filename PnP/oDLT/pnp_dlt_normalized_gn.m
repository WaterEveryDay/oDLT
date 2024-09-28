function [R, t] = pnp_dlt_normalized_gn(X, U, K)
% X = npoints x 3
% U = npoints x 3

% n = size(X, 1);
% n_small = min(n, min( max(30, floor(n/10)), 100) );
% [Rp, Rt] = pnp_dlt_normalized(X(1:n_small,:), U(1:n_small,:), K);

[R, t] = pnp_dlt_normalized(X, U, K);
[R, t] = optimize_pose_gn(X, U, K, R, t);
end
