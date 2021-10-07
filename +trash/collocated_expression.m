function obj = collocated_expression(N, x, y)
obj = yop.lagrange_signal.empty(1,0);
for n=1:N
    obj(n) = yop.lagrange_signal(x, y{n});
end
obj(N+1) = yop.lagrange_signal(0, y{N+1});
end