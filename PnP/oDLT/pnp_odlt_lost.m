function [R, t] = pnp_odlt_lost(X, U, K, Sig_uu)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with optimal DLT and refine the position with LOST
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

% solve optimal DLT for rotation
[R, t] = pnp_odlt(X, U, K, Sig_uu);

% solve optimal DLT for position
K_inv = invert_K(K);
r = refine_LOST(K_inv*U', X', R, K_inv*Sig_uu*K_inv', -R'*t);
if ~isnan(r)
    t = -R*r;
end
end