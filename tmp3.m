import yop.*

[t, t0, tf] = independent('t');
x = state('x', 2);
u = control('u');

s = x(1);
v = x(2);
a = u;

dx1 = der(s) == v;
