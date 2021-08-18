v = casadi.MX.sym('v');
b = casadi.MX.sym('b');

ff = casadi.Function('ff', {v, b}, {b});

ff(2, 1)