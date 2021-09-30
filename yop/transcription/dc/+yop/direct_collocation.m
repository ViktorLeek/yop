function nlp = direct_collocation(ocp, N, d, cp)
tau = [0, casadi.collocation_points(d, cp)];

t0 = yop.cx('t0');
tf = yop.cx('tf');
t = yop.collocated_time(t0, tf, N);
x = yop.collocated_state(n_x(ocp), N, tau);
z = [];
u = yop.collocated_control(n_u(ocp), N);
p = yop.cx('p', n_p(ocp));
dt = (tf - t0)/N;

[~, T0, Tf] = fixed_horizon(ocp);
[tps, ints] = yop.param_special_nodes(ocp.special_nodes, ocp.n_tp, ...
    ocp.n_int, N, tau, dt, T0, Tf, t, x, u, p);
ders = [];

J = ocp.objective.fn(p, tps, ints);

g = yop.discretize_ode(ocp, N, tau, dt, t, x, u, p);
for pc = ocp.equality_constraints
    disc = yop.parameterize_expression(pc,N,tau,t,x,u,p,tps,ints,ders);
    g = [g; disc(:)];
end

h = [];
for pc = ocp.inequality_constraints
    disc = yop.parameterize_expression(pc,N,tau,t,x,u,p,tps,ints,ders);
    h = [h; disc(:)];
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
nlp.z = z;
nlp.u = u;
nlp.p = p;
nlp.N = N;
nlp.d = d;
nlp.tau = tau;
nlp.cp = cp;
nlp.dt = dt;
end