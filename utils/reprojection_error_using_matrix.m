function err = reprojection_error_using_matrix(U, X, Phat)
% U = 3xnmeas where z component is 1
% X = 4xnmeas
% H = camera matrix
n_meas = size(X, 2);
Uhat = Phat*[X; ones(1, n_meas)];
Uhat = Uhat ./ Uhat(3,:);

eps = vecnorm(Uhat - U, 2, 1);
err = mean(eps);
end