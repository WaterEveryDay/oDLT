function ax = cross_mat( a )
%%%%%
%
% Function Name: cross_mat.m
%
%%%%%

ax = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
