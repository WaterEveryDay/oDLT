function [R, t] = pnp_epnp(X, U, K)
nmeas = size(X, 1);
[R, t] = efficient_pnp([X, ones(nmeas, 1)], U, K);
end
