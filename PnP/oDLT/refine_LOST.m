function [rI] = refine_LOST(xyb_all, pI_all, T_ItoC, R_xx, rI)
% Sebastien Henry and John Christian
% June 2022
% 
% The LOST algorithm is a weighted version of the DLT and it allows to
% perform global and statistically optimal triangulation without
% iteration. This code triangulates the position of the object assuming
% isotropic noise. Please cite Ref. [1] if using this algorithm.
%
% References
% - [1] S. Henry and J. A. Christian. Absolute Triangulation Algorithms for
%   Space Exploration. JGCD (2023). https://doi.org/10.2514/1.G006989
%
% Inputs:
% - xyb_all (3xn_meas): homogeneous measurements in the camera
%   frame, normalized such that x(3) = 1
% - pI_all (3xn_meas): 
%   positions of the observed objects in the inertial frame (resection)
%   OR
%   positions of of the cameras (intersection)
% - T_ItoC_all (3x3xn_meas): Rotation matrices from inertial to camera
%   frame for each of the respective measurements
% - R_xx_all (3x3xn_meas): Covariance matrices of the respective
%   measurements. ! this script assumes isotropic noise and only cares
%   about the (1, 1, ii) entry of R_xx = sig_xi^2.
%
% Outputs:
% - rI: position of the triangulated point in the inertial frame.
%   observer's position (resection)
%   OR 
%   observed position (intersection)
% - Prr: covariance of the inertial position 

% obtain the number of measurements
[~, n_meas] = size(xyb_all); 

% preallocation to speed-up code
H = zeros(2*n_meas, 3);
y = zeros(2*n_meas, 1);
sigx = sqrt(R_xx(1, 1)); % ! this code only works for isotropic noise

% extract matrix rotation 
t11 = T_ItoC(1, 1);
t12 = T_ItoC(1, 2);
t13 = T_ItoC(1, 3);

t21 = T_ItoC(2, 1);
t22 = T_ItoC(2, 2);
t23 = T_ItoC(2, 3);

t31 = T_ItoC(3, 1);
t32 = T_ItoC(3, 2);
t33 = T_ItoC(3, 3);

% extract image plane measurements
u_all = xyb_all(1,:)';
v_all = xyb_all(2,:)';

% compute LOST weights
weights = (vecnorm(xyb_all, 2, 1) ./ (sigx * vecnorm(pI_all - rI, 2, 1)))';

% compute DLT system with applied weights
% - add the rows with respect to measurements v
H(1:n_meas, :)          = weights.*[-t21+v_all*t31, -t22+v_all*t32, -t23+v_all*t33]; 
% - add the rows with respect to measurements u
H(n_meas+1:2*n_meas, :) = weights.*[ t11-u_all*t31,  t12-u_all*t32,  t13-u_all*t33];

% compute the y vector of the DLT system
y(1:n_meas) = sum(H(1:n_meas, :) .* pI_all', 2);
y(n_meas+1:2*n_meas) = sum(H(n_meas+1:2*n_meas, :) .* pI_all', 2);

% solve the system
rI = H \ y;
end
