function [R, t] = pnp_dlt(X, U, K)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with the regular DLT
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - X (nx3): 3D coordinates of the points in the world frame
% - U (nx3): Pixel coordinates of the corresponding points
% - K (3x3): Camera intrinsic calibration matrix
%
% Outputs:
% - R (3x3): Rotation matrix from world to camera frame
% - t (3x1): Translation vector representing the camera position in the world frame

% build the system
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
