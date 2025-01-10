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

    W_vec = sqrt(W(:));

    % Convert the initial guess matrix A into an angle vector for optimization
    [U, ~, V] = svd(A);
    R0 = U * diag([1, 1, det(U * V')]) * V';

    % Compute linear deviation for optimal procrustes problem
    % Linearize dcm = dcm0 * (eye - crossmat(dphi))
    eps = R0(:) - A(:);

    Jphi = - kron(eye(3), R0) * [...
                    0, 0, 0;
                    0, 0, 1;
                    0, -1, 0;
                    0, 0, -1;
                    0, 0, 0;
                    1, 0, 0;
                    0, 1, 0;
                    -1, 0, 0;
                    0, 0, 0;
                    ];

    dphi = -(W_vec.*Jphi) \ (W_vec.*eps);

    R = R0 * rotVecToMat(dphi);
end