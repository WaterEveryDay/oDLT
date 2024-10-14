addpath ../PnP/oDLT/
theta = 0*pi/180;
gamma = 0*pi/180;
R_I2C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)] ...
    * [cos(gamma), sin(gamma), 0; -sin(gamma), cos(gamma), 0; 0, 0, 1];

% camera position
r = [0; 0; 0];

% camera intrinsics
dx = 0.001250;
dy = dx;
K = [1/dx, 0, 320; 0, 1/dy, 240; 0, 0, 1]; %eye(3);

verbose = true;

mode = "uncentered";
precise_timing = true;

% default parameters
n_meas = 50;
sig_u = 1;
outlier_percentage = 0;

% parameter changing, select
var_changing = "n";

% Monte Carlo
nmc = 10000;

switch mode
    case "centered"
        bbox = [-2, -2, 4; 2, 2, 8];
    case "uncentered"
        bbox = [1, 1, 4; 2, 2, 8];
    case "own1"
        %bbox = [-5, -5, 10; 5, 5, 15];
        bbox = [-5, -5, 10; 5, 5, 15];
end

%% define pixel boundaries
ubox = [0, 0; 640, 480];

Sigma_uu = diag([sig_u^2, sig_u^2, 0]);

H = K*R_I2C*[eye(3), -r];

% generate 3D points randomly
X = unifrnd(repmat(bbox(1,:)', 1, n_meas), ...
    repmat(bbox(2,:)',1, n_meas), 3, n_meas);

% project 3D points to pixel
U = H*[X; ones(1, n_meas)];
U = U./U(3,:);

hvec = zeros(nmc, 12);
% iterate over the monte carlo
for i = 1:nmc


    % make pixel measurement noisy
    Utilde = U + mvnrnd([0; 0; 0], Sigma_uu, n_meas)';
    [h, cov_h] = pnp_odlt_matrix(X', Utilde', Sigma_uu);
    hvec(i,:) = h';

end
[h, cov_h] = pnp_odlt_matrix(X', U', Sigma_uu);
cov_closed_form = cov_h ./ cov_h(1, 1)
cov_mc = cov(hvec);
cov_mc = cov_mc ./ cov_mc(1, 1)

function [h, cov_h] = ... %hbar_uncalibrated, cov_hbar_uncalibrated] = ...
    pnp_odlt_matrix(X, U, Sig_uu)
% X = 3xn: 3D points in world
% U = 3xn: 2D pixel measurements
% K = 3x3: Calibration matrix
% Sig_uu = 3x3: Covariance of pixel measurements

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

n = size(X, 1);

% normalize the measurements and 3D points for a better behaviour
[X_scaled, T_X, ~] = normalize_points3D(X);
[U_scaled, T_U, T_U_inv] = normalize_meas(U);

% build the matrix for DLT system
A = build_dlt_system(X_scaled, U_scaled);

% solve for an initial guess using a subset of the system
n_small = min(2*n, min( max(30, floor(n/10)), 100) );

% rng(1)
[~,~,V] = svd(A(randperm(2*n, n_small), :), "econ");

h = V(:,end);
H_init = reshape(h, 3, 4);

% compute weights for the DLT system
sig_u = T_U(1,1)*sqrt(Sig_uu(1,1));
X_proj = H_init * [X_scaled, ones(n, 1)]';
scale_3 = X_proj(3,:);
inv_scale_norm_sig = 1./ (scale_3 * sig_u);

% apply weights to the DLT system
A = kron(inv_scale_norm_sig', [1; 1]) .* A;

% solve the system by finding smallest eigen vector of A^T A
[~, Sig, V] = svd(A, 'econ');
h = V(:,end);
% cov_h = inv(A'*A);
d = diag(Sig);

if d(end) / d(1) < 1e-7
    d(end) = 1e-7*d(1);
end

cov_h = V * diag(1./d.^2) * V';
end