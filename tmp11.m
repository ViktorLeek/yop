
d = 2;
tau = casadi.collocation_points(d, 'legendre');
[C,D,B] = casadi.collocation_coeff(tau);
C = full(C)
D = full(D);
B = full(B);

%%
L = yop.lagrange_polynomial([0 tau], sym('tau', [1, d+1]));
% L.l
L.eval_basis(1)

Ld = L.differentiate();
Ld.eval_basis([0 tau]);

Li = L.integrate();
Li.eval_basis([0 tau]);

%%
tmp = yop.lagrange_basis(0);
l.eval_basis(0)