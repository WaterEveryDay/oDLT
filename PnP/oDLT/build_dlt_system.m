function A = build_dlt_system(X, U)
    % X = 3 x n: 3D points
    % U = 2 x n or 3 x n: pixel measurements
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
    A(1:2:2*n,:) = [x,  zero, -ux, y,  zero, -uy, z, zero, -uz, one, zero, -u]; 
    A(2:2:2*n,:) = [zero, -x, vx, zero, -y, vy, zero, -z, vz, zero, -one, v];
end

