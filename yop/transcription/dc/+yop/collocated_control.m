function obj = collocated_control(nu ,N)
y = cell(N+1,1);
for n=1:N
    y{n} = casadi.MX.sym(['u_' num2str(n)], nu, 1);
end
y{N+1} = y{N};
obj = yop.collocated_expression(N, 0, y);

% Bad code! But avoids other functions having insight into the structure of
% the lagrange polynomials.
for ok = obj
    ok.exclude_last = true;
end
end