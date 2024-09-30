function [R, t, K] = decompose_P(H)

H_inv = inv(H(1:3,1:3));

Rot = [-1, 0, 0; 0, -1, 0; 0, 0, 1];
% recover the inverse rotation and inverse of K
[R_transpose, K_inv] = qr(H_inv);

% recover the calibration matrix
K = inv(K_inv);

% scale such that the last element of the calibration matrix is 1
scale = K(3,3);
K = K*Rot/scale;

% recover extrinsics
R = Rot*R_transpose';
r = - H_inv * H(:,4);
t = - R * r;
end