function nlp = direct_collocation(ocp, N, d, cp)
tau = [0, casadi.collocation_points(d, cp)];

[is_fixed, t0, tf] = fixed_horizon(ocp);
if ~is_fixed
    t0 = casadi.MX.sym('t0');
    tf = casadi.MX.sym('tf');
end
dt = (tf - t0)/N;

t = yop.collocated_time(t0, tf, N);
x = yop.collocated_state(n_x(ocp), N, tau);
u = yop.collocated_control(n_u(ocp), N);
p = casadi.MX.sym('p', n_p(ocp));

[tps, ints] = yop.param_special_nodes(ocp.special_nodes, ocp.n_tp, ...
    ocp.n_int, N, tau, dt, t0, tf, t, x, u, p);

J = ocp.objective.fn(tps, ints);

g = yop.discretize_ode(ocp, N, tau, dt, t, x, u, p);
for pc = ocp.equality_constraints
    g = [g; yop.parameterize_expression(pc, N, tau, t, x, u, p, tps, ints)];
end

h = [];
for pc = ocp.inequality_constraints
    h = [h; yop.parameterize_expression(pc, N, tau, t, x, u, p, tps, ints)];
end

[w_lb, w_ub] = yop.box_bnd(N, d, ocp);

nlp = struct;
nlp.J = J;
nlp.w = vertcat(t0, tf, vec(x), vec(u), p);
nlp.w_ub = w_ub;
nlp.w_lb = w_lb;
nlp.g = g;
nlp.g_ub = zeros(size(g));
nlp.g_lb = zeros(size(g));
nlp.h = h;
nlp.h_ub = zeros(size(h));
nlp.h_lb = -inf(size(h));
nlp.t = t;
nlp.x = x;
nlp.u = u;
nlp.p = p;
nlp.N = N;
nlp.d = d;
nlp.tau = tau;
nlp.cp = cp;
nlp.dt = dt;
end