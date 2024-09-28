function [R, t] = pnp_opnp(X, U, K)
[R, t] = OPnP(X', K\U', 'polish');
end
