import yop.*

[t, t0, tf] = independent('t');
x = state('x', 2);
u = control('u');

p = x(1);
v = x(2);
a = u;

ocp = yop.ocp();
ocp.min( 0.5 * int(a^2) );
ocp.st(...
    t0 == 0, tf == 1, ...
    der(p) == v, ...
    der(v) == a, ...
    p(t0) == p(tf) == 0, ...
    v(t0) == -v(tf) == 1, ...
    p <= 1/9 ...
    );

ocp.parse_variables();
ocp.parse_constraints();