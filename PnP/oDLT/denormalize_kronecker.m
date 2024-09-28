function [h, cov_h] = denormalize_kronecker(h, cov_h, K_inv, T_U_inv, T_X)
    % Step 1: Create the transformation matrix
    M = kron(T_X', K_inv * T_U_inv);
    
    % Step 2: Apply the transformation to h
    h = M * h;
    
    % Step 3: Transform the covariance matrix
    cov_h = M * cov_h * M';
end

