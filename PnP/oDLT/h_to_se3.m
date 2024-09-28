function [R, t] = h_to_se3(h)
    % inputs
    % h: 12x1: wrongly scaled and unorthogonalized camera matrix H = [R, -R'*r]

    % output: 
    
    % extract the rotation matrix and compute the scaling
    R_init = reshape(h(1:9), 3, 3);
    detR = det(R_init);
    scale = sign(detR) * abs(detR)^(1/3);
    
    % scale h
    h = h/scale;
    R_init = reshape(h(1:9), 3, 3);
    [U, ~, V] = svd(R_init);
    
    % retrieve SO3 matrix and translation
    R = U *diag([1, 1, det(U*V')])* V';
    t = h(10:12);
end