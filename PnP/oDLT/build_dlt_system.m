function A = build_dlt_system(X, U)
% Sebastien Henry
% Efficiently build the DLT system, accounting for the third redundant row.
%
% References:
% - [1] S. Henry and J. A. Christian. Optimal DLT-based Solutions for the Perspective-n-Point. (2024).
%
% Inputs:
% X = nx3: 3D points in the world coordinates
% U = nx2 or nx3: 2D pixel measurements (non-homogeneous or homogeneous with z = 1)
%
% Outputs:
% A = 3nx3: Matrix incorporating the DLT constraint

n = size(X, 1);

% extract columns
x = X(:,1);
y = X(:,2);
z = X(:,3);
u = U(:,1);
v = U(:,2);

% pre-compute multiplications
ux = u .* x;
uy = u .* y;
uz = u .* z;
vx = v .* x;
vy = v .* y;
vz = v .* z;

% pre-compute matrices of zeros and ones
zero = zeros(n, 1);
one = ones(n, 1);

% fill A using vectorized operations
A = zeros(2*n, 12);
A(1:2:2*n,:) = [zero, -x, vx, zero, -y, vy, zero, -z, vz, zero, -one, v];
A(2:2:2*n,:) = [x,  zero, -ux, y,  zero, -uy, z, zero, -uz, one, zero, -u];
end

