function [R, t] = h_to_se3_with_cov(h, cov_h)
    % inputs
    % h: 12x1: wrongly scaled and unorthogonalized camera matrix H = [R, -R'*r]
    % cov_h: 12x12: covariance of h

    % output: 
    
    % extract the rotation matrix and compute the scaling
    R_init = reshape(h(1:9), 3, 3);
    detR = det(R_init);
    scale = sign(detR) * abs(detR)^(1/3);
    
    % scale h
    h = h/scale;
    R_init = reshape(h(1:9), 3, 3);
    
    % check whether cov_h is suitable 
    if rank(cov_h) < 12
        W = eye(3);
    else
        W = scale^2 * (cov_h)^(-1);
        W = reshape(diag(W(1:9, 1:9)), 3, 3);
        W = diag(abs(sum(W, 1)));
    end
    [U, ~, V] = svd(R_init*W);
    R = U *diag([1, 1, det(U*V')])* V';    
    t = h(10:12);
end