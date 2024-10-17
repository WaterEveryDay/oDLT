syms u1 u2 p1 p2 p3 P11 P12 P13 P14 P21 P22 P23 P24 P31 P32 P33 P34 sigu
assume([u1 u2 p1 p2 p3 P11 P12 P13 P14 P21 P22 P23 P24 P31 P32 P33 P34 sigu], 'real')
ubar = [u1; u2; 1];
pbar = [p1; p2; p3; 1];
k = [1; 0; 0];
S = [1 0 0; 0 1 0];
P = [P11 P12 P13 P14; P21 P22 P23 P24; P31 P32 P33 P34];

Bi = S * cross_mat(ubar)^2 / (sigu * norm(ubar)^2 * k' * P * pbar);
Ai = kron(pbar', cross_mat(ubar));



qi = 1/(sigu * k' * P * pbar);

exp1 = Bi * Ai;
exp2 = -qi * S * Ai;
simplify(exp1 - exp2)