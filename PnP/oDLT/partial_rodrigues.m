function dRdg = partial_rodrigues(g)
% Sebastien Henry

% Rodrigues parameters are such that tan(theta/2)*e, where e is euler
% vector
% compute the jacobian between rodrigues and rotation matrix

% retrieve the vectorized rotation matrix
g1 = g(1); g2 = g(2); g3 = g(3);
denom = 1+(g1^2+g2^2+g3^2);

Rvec = 1/denom * ...
       [1 + g1^2 - g2^2 - g3^2;
        2*(g2*g1 - g3);
        2*(g3*g1 + g2);
        2*(g1*g2 + g3);
        1 - g1^2 + g2^2 - g3^2;
        2*(g3*g2 - g1);
        2*(g1*g3 - g2);
        2*(g2*g3 + g1);
        1 - g1^2 - g2^2 + g3^2];

% retrieve derivative 
dRvecdg1 = 2/denom * ...
          [ g1, -g2, -g3;
            g2,  g1,  -1;
            g3,   1,  g1;
            g2, g1,   1;
           -g1, g2, -g3;
            -1, g3,  g2;
            g3,  -1, g1;
             1,  g3, g2;
           -g1, -g2, g3];

dRvecdg2 = -(2/denom) * Rvec*g';


dRdg = dRvecdg1 + dRvecdg2;





% num = ((1 -(g1^2+g2^2+g3^2))*eye(3) + 2* (g*g') - 2*cross_mat(g));
% denom = 1+(g1^2+g2^2+g3^2);
% 
% dnumddg_times_denom = (2/denom) * ...
% [g1, -g2, -g3;
% g2,  g1,  -1;
% g3,   1,  g1;
%  g2, g1,   1;
% -g1, g2, -g3;
%  -1, g3,  g2;
%  g3,  -1, g1;
%   1,  g3, g2;
% -g1, -g2, g3];
% 
% ddenomdg_times_num = ...
% -(2/denom^2)*[num(:,1)*g'; num(:,2)*g'; num(:,3)*g'];
% 
% dRdg2 = dnumddg_times_denom + ddenomdg_times_num;

end

