function R = rodrigues2mat(g)
% cross_g = cross_mat(g);
% R = eye(3) + 2 * (cross_g*cross_g - cross_g) / (1 + norm(g)^2) ;

g1 = g(1); g2 = g(2); g3 = g(3);
denom = 1+(g1^2+g2^2+g3^2);

R = 1/denom * ...
    [1 + g1^2 - g2^2 - g3^2, 2*(g1*g2 + g3)        , 2*(g1*g3 - g2);
     2*(g2*g1 - g3)        , 1 - g1^2 + g2^2 - g3^2, 2*(g2*g3 + g1);
     2*(g3*g1 + g2)        , 2*(g3*g2 - g1)        , 1 - g1^2 - g2^2 + g3^2];
end

