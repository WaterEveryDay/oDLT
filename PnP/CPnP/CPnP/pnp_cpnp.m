function [R, t] = pnp_cpnp_gn(X, U, K)

% format data
fx = K(1,1);
fy = K(2,2);
u0 = K(1,3);
v0 = K(2,3);
U = U(:,1:2)';

% run CPnP
[R, t] = CPnP(X', U,fx,fy,u0,v0);

end

