function R = weighted_orthogonal_procrustes_gn(A, W)
    % Sebastien Henry
    % 
    % Function to find the closest orthogonal matrix B to A in SO(3)
    % with element-wise weights W using lsqnonlin
    
    % Input:
    % A - 3x3 input matrix
    % W - 3x3 weight matrix
    % Output:
    % B_opt - Optimal orthogonal matrix B closest to A in SO(3)
    
    % Ensure that A and W are 3x3 matrices
    assert(isequal(size(A), [3, 3]), 'A must be a 3x3 matrix.');
    assert(isequal(size(W), [3, 3]), 'W must be a 3x3 matrix.');
    
    A_vec = A(:);
    W_vec = sqrt(W(:));

    % Convert the initial guess matrix A into an angle vector for optimization
    [U, ~, V] = svd(A);
    R = U * diag([1, 1, det(U * V')]) * V';

    g = mat2rodrigues(R);

    iter = 0;
    max_iter = 1;

    while iter < max_iter
        % compute partial
        dRdg = partial_rodrigues(g);

        % compute residual between both vectors
        epsilon = A_vec - R(:);

        % update
        dg = (W_vec.*dRdg) \ (W_vec.*epsilon);
        
        g = g + dg;
        R = rodrigues2mat(g);
        iter = iter + 1;
    end
end
