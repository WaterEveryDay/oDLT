function dtheta = err_DCM(R, Rhat)
% Crassidis & Markley p 224
dtheta = 2*asin( norm(Rhat-R, 'fro')/sqrt(8) );
end
