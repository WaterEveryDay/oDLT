function [R, t] = pnp_dlt_normalized_gn(X, U, K)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Solve the PnP with normalized DLT and then refine with GN
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


[R, t] = pnp_dlt_normalized(X, U, K);
[R, t] = optimize_pose_gn(X, U, K, R, t);
end
