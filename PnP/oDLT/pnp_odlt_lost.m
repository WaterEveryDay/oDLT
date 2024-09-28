function [R, t] = pnp_opt_dlt_lost(X, U, K, Sig_uu)
% X (nx3) 3D points in world
% U (nx3) 2D pixel measurements
% K (3x3) Calibration matrix
% Sig_uu (3x3) Pixel covariance (optional)

% arguments
%     X(:,3) double
%     U(:,3) double
%     K(3,3) double
%     Sig_uu(3,3) double = [1, 0, 0; 0, 1, 0; 0, 0, 0];
% end
if nargin < 4
    Sig_uu = [1, 0, 0; 0, 1, 0; 0, 0, 0];
end

% solve optimal DLT for rotation
[R, t] = pnp_odlt(X, U, K, Sig_uu);

% solve optimal DLT for position
K_inv = invert_K(K);
r = refine_LOST(K_inv*U', X', R, K_inv*Sig_uu*K_inv', -R'*t);
if ~isnan(r)
    t = -R*r;
end
end