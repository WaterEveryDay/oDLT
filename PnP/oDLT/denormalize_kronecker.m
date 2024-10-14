function [h, inv_cov_h] = denormalize_kronecker_inv(h, inv_cov_h, K_inv, T_U_inv, T_X)
    % Step 1: Create the transformation matrix
    M = kron(T_X', K_inv * T_U_inv);
    
    % Step 2: Apply the transformation to h
    h = M * h;
    
    M_inv = inv(M);
    % Step 3: Transform the covariance matrix
    inv_cov_h = M_inv' * inv_cov_h * M_inv;
end

