function obj = collocated_time(t0, tf, N)
x = [0, 1];
dt = (tf-t0)/N;
y = cell(N+1,1);
for n=1:N
    y{n} = t0 + dt*[n-1, n];
end
y{N+1} = tf;
obj = yop.collocated_expression(N, x, y);
end

