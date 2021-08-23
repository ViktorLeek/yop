import yop.*

[t, t0, tf] = independent('t');
x = state('x', 1, 2);
u = control('u');
z = algebraic('z', 1, 4);

tmp = [x(1), x(2), x, u];
tmp2 = [tmp(4), tmp(3), tmp(1), tmp(2), tmp(4)];
[b, v] = isa_variable([tmp2(2:end), z(2:3)+2])
%%
tmp = [true, false, true, true, false];
all([true, true, false])

x = yop.state('x', 3, 3);
tmp = x(1:2, 1:2);

sz = size(tmp.node);
tmp.node.m_value = reshape(1:prod(sz), sz);
forward(tmp)


%%
x = yop.state('x', 3);
y = yop.state('y', 3);

x1 = x(1);
x2 = x(2);
y3 = y(3);


m = [x1, x1+x2, x1; x2+1, x2, x2; y3, der(y3), y3];
v = [x; x2; y3-2; x1+1];

v12 = v(end-1:end);

v12 <= 2;

isa_variable([x1, x2])