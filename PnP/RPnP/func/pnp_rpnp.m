function [R, t] = pnp_rpnp(X, U, K)

S = [1, 0, 0; 0, 1, 0];
U = S*invert_K(K)*U';
[R, t] = RPnP(X', U);
end

