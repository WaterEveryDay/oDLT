function [R, t] = pnp_epnp_gn(X, U, K)
nmeas = size(X, 1);
[R, t, ~, ~] = efficient_pnp_gauss([X, ones(nmeas, 1)],U,K);
end
