function obj = collocated_state(nx ,N, tau)
y = cell(N+1,1);
for n=1:N
    y{n} = casadi.MX.sym(['x_' num2str(n)], nx, length(tau));
end
y{N+1} = casadi.MX.sym(['x_' num2str(N+1)], nx);
obj = yop.collocated_expression(N, tau, y);

end

