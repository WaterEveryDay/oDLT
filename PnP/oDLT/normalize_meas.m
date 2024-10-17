function [Utilde, T, T_inv] = normalize_meas(U)
% Authors: Sebastien Henry
% Last Modified: October 2024
%
% Normalize the pixel measurement and compute the similarity transformation matrices
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% - U (nx3): Pixel measurements [u_i, v_i, 1]
%
% Outputs:
% - Utilde (nx3): points that are resenctered and scaled
% - T (3x3): similarity transformation matrix
% - T_inv (3x3): inverse similarity transformation matrix

% Step 1: Compute the centroid of the points
centroid = mean(U, 1);  % [c_x, c_y]

% Step 2: Translate the points so the centroid is at the origin
U_translated = U - centroid;  % Subtract centroid from each point

% Step 3: Compute the average distance of the translated points from the origin
d2_avg = mean(sum(U_translated(:,1:2).^2, 2));

% Step 4: Scale the points so that the average distance is sqrt(2)
scale_factor = sqrt(2 / d2_avg);
Utilde = scale_factor * U_translated;
Utilde(:,3) = ones(size(U, 1), 1);

% Step 5: Construct the similarity transformation matrix T
% T includes translation and scaling
T = [scale_factor, 0, -scale_factor * centroid(1);
    0, scale_factor, -scale_factor * centroid(2);
    0, 0, 1];

T_inv = [1/scale_factor, 0, centroid(1);
    0, 1/scale_factor, centroid(2);
    0, 0, 1];

end