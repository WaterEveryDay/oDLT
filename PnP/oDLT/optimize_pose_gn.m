% CPnP: a consistent PnP solver
% Inputs: s - a 3*n matrix whose i-th column is the coordinates (in the world frame) of the i-th 3D point
%        Psens_2D - a 2*n matrix whose i-th column is the coordinates of the 2D projection of the i-th 3D point
%        fx, fy, u0, v0 - intrinsics of the camera, corresponding to the intrinsic matrix K=[fx 0 u0;0 fy v0;0 0 1]
%
%Outputs: R - the estimate of the rotation matrix in the first step
%         t - the estimate of the translation vector in the first step
%         R_GN - the refined estimate of the rotation matrix with Gauss-Newton iterations
%         t_GN - the refined estimate of the translation vector with Gauss-Newton iterations

% Copyright <2022>  <Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, Junfeng Wu>
% Guangyang Zeng, SLAMLab-CUHKSZ, September 2022
% zengguangyang@cuhk.edu.cn, https://github.com/SLAMLab-CUHKSZ 
% paper information: Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, and Junfeng Wu. 
%                    CPnP: Consistent Pose Estimator for Perspective-n-Point Problem with Bias Elimination, 
%                    IEEE International Conference on Robotics and Automation (ICRA), London, UK, May 2023.

% function [R,t,R_GN,t_GN]=CPnP(s,Psens_2D,fx,fy,u0,v0)
function [R, t] = optimize_pose_gn(X, U, K, R, t)
% Modified October 2024: Sebastien Henry - Gauss Newton adapted from CPnP paper

N = size(X, 1);
S = [1 0 0;0 1 0];
W = diag([K(1,1), K(2,2)]);
WS=W*S;

Psens_2D = U(:,1:2)' - [K(1,3); K(2,3)];
obs=Psens_2D(:);
s = X';

e3 = [0;0;1];

iter = 0;
max_iter = 3;
while iter < max_iter
   
    J = zeros(2*N, 6);
    g=WS* (R*s+t);
    h=e3'* (R*s+t);
    f=g./h;
    f=f(:);
    I3=eye(3);
    for k = 1:N
        J(2*k-1:2*k,:) = (((WS * h(k) - g(:,k)* e3') * [s(2,k)*R(:,3)-s(3,k)*R(:,2) s(3,k)*R(:,1)-s(1,k)*R(:,3) s(1,k)*R(:,2)-s(2,k)*R(:,1) I3]) )/ h(k)^2;
    end
    initial = [0;0;0;t];
    dpose = (J' * J)\ (J') * (obs - f);
    results = initial + dpose;
    X_GN = results(1:3);
    t = results(4:6);
    Xhat = [0 -X_GN(3) X_GN(2); X_GN(3) 0 -X_GN(1); -X_GN(2) X_GN(1) 0];
    R = R * expm(Xhat);
    iter = iter+1;
end


end