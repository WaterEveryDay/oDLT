function dRdg = partial_rodrigues(g)
% Sebastien Henry
% see attached test script for confirmation

g1 = g(1); g2 = g(2); g3 = g(3);

num = ((1 -(g1^2+g2^2+g3^2))*eye(3) + 2* (g*g') - 2*cross_mat(g));
denom = 1+(g1^2+g2^2+g3^2);

dnumddg_times_denom = (2/denom) * ...
[g1, -g2, -g3;
g2,  g1,  -1;
g3,   1,  g1;
 g2, g1,   1;
-g1, g2, -g3;
 -1, g3,  g2;
 g3,  -1, g1;
  1,  g3, g2;
-g1, -g2, g3];

ddenomdg_times_num = ...
-(2/denom^2)*[num(:,1)*g'; num(:,2)*g'; num(:,3)*g'];

dRdg = dnumddg_times_denom + ddenomdg_times_num;

end

