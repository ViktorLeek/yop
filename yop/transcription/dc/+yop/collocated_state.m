function obj = collocated_state(name, sig_dim ,N, polydeg, points)

x = [0, casadi.collocation_points(polydeg, points)];

obj = yop.lagrange_polynomial.empty(1,0);
for n=1:N
    y = casadi.MX.sym(label(name,n), sig_dim, polydeg+1);
    obj(n) = yop.lagrange_polynomial(x, y);
end
y = casadi.MX.sym(label(name, N+1), sig_dim);
obj(N+1) = yop.lagrange_polynomial(0, y);

end

function l = label(name, n)
l = [name, '_' num2str(n)];
end
