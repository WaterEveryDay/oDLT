function ax = cross_mat( a )
%%%%%
%
% Function Name: cross_mat.m
%
%%%%%

ax = zeros(3,3);

ax(3,2) =  a(1);
ax(2,3) = -a(1);

ax(1,3) =  a(2);
ax(3,1) = -a(2);

ax(2,1) =  a(3);
ax(1,2) = -a(3);
