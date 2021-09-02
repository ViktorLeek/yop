function rk4 = rk4_dms_expr(ode, q, nx, nu, np, nw, steps)

if nargin == 5
    steps = yop.defaults.rk4_steps;
end

p0  = casadi.MX.sym('t0'); % Parameter t0 for entire problem
pf  = casadi.MX.sym('tf'); % Parameter tf for entire problem
t0  = casadi.MX.sym('t0');
tf  = casadi.MX.sym('tf');
tt  = casadi.MX.sym('t');
x0 = casadi.MX.sym('x0', nx);
u  = casadi.MX.sym('U', nu);
p  = casadi.MX.sym('P', np);
w = casadi.MX.sym('w', nw);
t = t0;
x = x0;

F = casadi.Function('F', {p0, pf, tt, x, u, p, w}, ...
    {ode(p0, pf, tt, x, u, p), q(p0, pf, tt, x, u, p, w)});

h = (tf-t0)/steps;
Q = 0;
for j=1:steps
    [k1, q1] = F(p0, pf, t      , x           , u, p, w);
    [k2, q2] = F(p0, pf, t + h/2, x + h/2 * k1, u, p, w);
    [k3, q3] = F(p0, pf, t + h/2, x + h/2 * k2, u, p, w);
    [k4, q4] = F(p0, pf, t + h  , x + h   * k3, u, p, w);
    t = t + h;
    x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    Q = Q + h/6*(q1 + 2*q2 + 2*q3 + q4);
end
rk4 = casadi.Function('Q', {t0, tf, p0, pf, x0, u, p, w}, {Q});

end