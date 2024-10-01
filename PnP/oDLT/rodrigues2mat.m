function R = rodrigues2mat(g)
   cross_g = cross_mat(g);
   R = eye(3) + 2 * (cross_g*cross_g - cross_g) / (1 + norm(g)^2) ;
end

