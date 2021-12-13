function nlp = direct_collocation2(ocp, N, d, cp)

% If horizon is not fixed, T0 and Tf are nonsense. But since they are only
% to be used for a fixed horizon, that is OK.
[~, T0, Tf] = ocp.fixed_horizon();

args = arguments(ocp, N, d, cp, T0, Tf);
args = param_special_nodes(ocp, N, T0, Tf, args);

J = ocp.objective.fn(args.t0, args.tf, args.p, args.tps, args.ints);

g = discretize_dynamics(ocp, N, args);
g_ub = zeros(size(g));
g_lb = zeros(size(g));

g = [g; pointcon(ocp.point, args)];
g_ub = [g_ub; ocp.point.ub];
g_lb = [g_lb; ocp.point.lb];

g = [g; pathcon(ocp.path, N, args)];
g_ub = [g_ub; ocp.path.ub];
g_lb = [g_lb; ocp.path.lb];

g = [g; hard_pathcon(ocp.path_hard, N, args)];
g_ub = [g_ub; ocp.path_hard.ub];
g_lb = [g_lb; ocp.path_hard.lb];

for pk = ocp.path_ival
    gk = ival_pathcon(pk, N, args);
    g = [g; gk];
    g_ub = [g_ub; ocp.path_ival.ub*ones(size(gk))];
    g_lb = [g_lb; ocp.path_ival.lb*ones(size(gk))];
end

[w_lb, w_ub] = box_bnd(T0, Tf, N, args.tau, ocp);

w = vertcat( ...
    args.t0, ...
    args.tf, ...
    args.x.vec(), ...
    args.z.vec(), ...
    args.u.vec(), ...
    args.p.vec() ...
    );

t0 = casadi.Function('t0', {w}, {args.t0});
tf = casadi.Function('tf', {w}, {args.tf});
t  = casadi.Function('t' , {w}, {args.t.mat()});
x  = casadi.Function('x' , {w}, {args.x.mat()});
z  = casadi.Function('z' , {w}, {args.z.mat()});
u  = casadi.Function('u' , {w}, {args.u.mat()});
p  = casadi.Function('p' , {w}, {args.p});

nlp = struct;
nlp.J = J;
nlp.w = w;
nlp.w_ub = w_ub;
nlp.w_lb = w_lb;
nlp.g = g;
nlp.g_ub = g_ub;
nlp.g_lb = g_lb;
nlp.t0 = @(w) full(t0(w));
nlp.tf = @(w) full(tf(w));
nlp.t  = @(w) full(t(w));
nlp.x  = @(w) full(x(w));
nlp.z  = @(w) full(z(w));
nlp.u  = @(w) full(u(w));
nlp.p  = @(w) full(p(w));
end

function args = arguments(ocp, N, d, cp, T0, Tf)

t0 = yop.cx('t0');
tf = yop.cx('tf');
dt = (tf - t0)/N;
tau = full([0, casadi.collocation_points(d, cp)]);

t = yop.interpolating_poly.empty(N+1, 0);
t_n = t0;
for n=1:N
    t(n) = yop.interpolating_poly(tau, t_n + tau*dt, T0, Tf, N);
    t_n = t_n + dt;
end
t(N+1) = yop.interpolating_poly(0, tf, T0, Tf, N);

x = yop.interpolating_poly.empty(N+1,0);
for n=1:N
    x_n = yop.cx(['x_' num2str(n)], ocp.n_x, d+1);
    x(n) = yop.interpolating_poly(tau, x_n, T0, Tf, N);
end
x_n = yop.cx(['x_' num2str(N+1)], ocp.n_x);
x(N+1) = yop.interpolating_poly(0, x_n, T0, Tf, N);

z = yop.interpolating_poly.empty(N, 0);
for n=1:N
    z_n = yop.cx(['z_' num2str(n)], ocp.n_z, d);
    z(n) = yop.interpolating_poly(tau(2:end), z_n, T0, Tf, N);
end

u = yop.interpolating_poly.empty(N, 0);
for n=1:N
    u_n = yop.cx(['u_' num2str(n)], ocp.n_u);
    u(n) = yop.interpolating_poly(0, u_n, T0, Tf, N);
end

p = yop.cx('p', ocp.n_p);

args = struct;
args.t0 = t0;
args.tf = tf;
args.t = t;
args.x = x;
args.z = z;
args.u = u;
args.p = p;
args.dt = dt;
args.tau = tau;

end


function args = param_special_nodes(ocp, N, T0, Tf, args)
    % (ocp, N, tau, dt, t0, tf, t, x, z, u, p, T0, Tf)

tps = [];
ints = [];
ders = yop.interpolating_poly.empty(N, 0);
for n=1:N
    ders(n) = yop.interpolating_poly(args.tau, [], T0, Tf, N);
end

for node = ocp.snodes
    tmp_tp  = [tps;  zeros(ocp.n_tp  - length(tps), 1)];
    tmp_int = [ints; zeros(ocp.n_int - length(ints), 1)];
    
    args.tps = tmp_tp;
    args.ints = tmp_int;
    args.ders = ders;
    
    switch node.type
        case yop.ocp_expr.tp
            tp = parameterize_timepoint(node, args, ocp.n_der);
            tps = [tps; tp(:)];
            
        case yop.ocp_expr.int
            int = parameterize_integral(node, N, args, ocp.n_der);
            ints = [ints; int(:)];
            
        case yop.ocp_expr.der
            ders = parameterize_derivative(node, N, args, ocp.n_der);
    end
end

args.tps = tps;
args.ints = ints;
args.ders = ders;

end

function val = parameterize_timepoint(tp, args, n_der)
t0 = args.t0;
tf = args.tf;
tt = args.t.value(tp.timepoint);
xx = args.x.value(tp.timepoint);
zz = args.z.value(tp.timepoint);
uu = args.u.value(tp.timepoint); 
pp = args.p;
tps = args.tps;
int = args.ints;
dd = args.ders.value(tp.timepoint);
dd = [dd;  zeros(n_der  - length(dd), 1)];
val = tp.fn(t0, tf, tt, xx, zz, uu, pp, tps, int, dd);
end

function I = parameterize_integral(i, N, args, n_der)
I = 0;
t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;
for n=1:N
    yval = [];
    for r=1:length(args.tau)
        tt = args.t(n).y(r);
        xx = args.x(n).y(:, r);
        zz = args.z(n).evaluate(args.tau(r));
        uu = args.u(n).y;
        dd = args.ders(n).evaluate(args.tau(r));
        dd = [dd;  zeros(n_der  - length(dd), 1)];
        val_r = i.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd);
        yval = [yval, val_r(:)];
    end
    lp = yop.lagrange_polynomial(args.tau, yval).integrate();
    I = I + lp.evaluate(1)*args.dt;
end
end

function ders = parameterize_derivative(der, N, args, n_der)
t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;
for n=1:N
    yval = [];
    for r=1:length(args.tau)
        tt = args.t(n).y(r);
        xx = args.x(n).y(:, r);
        zz = args.z(n).evaluate(args.tau(r));
        uu = args.u(n).y;
        dd = args.ders(n).evaluate(args.tau(r));
        dd = [dd;  zeros(n_der  - length(dd), 1)];
        val_r = der.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd);
        yval = [yval, val_r(:)];
    end
    lp = yop.lagrange_polynomial(args.tau, yval);
    yn = lp.differentiate().evaluate(args.tau)/args.dt;
    ders(n).y = [ders(n).y; yn];
end
end

function g = discretize_dynamics(ocp, N, args)
g = []; % Equality constraints from discretization
t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;
for n=1:N % Dynamics
    dx = args.x(n).differentiate();
    for r=2:length(args.tau)
        dxr = dx.evaluate(args.tau(r));
        tt = args.t(n).y(r);
        xx = args.x(n).y(:, r);
        zz = args.z(n).y(:, r-1); % only has d parameters
        uu = args.u(n).y(:);
        dd = args.ders(n).evaluate(args.tau(r));
        f = ocp.ode.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd);
        a = ocp.alg.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd);
        g = [g; (dxr - args.dt*f); a];
    end
end

for n=1:N % Continuity
    gn = args.x(n).evaluate(1) - args.x(n+1).evaluate(0);
    g = [g; gn];
end
end

function disc = pointcon(expr, args)
t0 = args.t0;
tf = args.tf;
tt = args.t(1).evaluate(0);
xx = args.x(1).evaluate(0);
zz = args.z(1).evaluate(0);
uu = args.u(1).evaluate(0);
pp = args.p;
tp = args.tps;
int = args.ints;
dd = args.ders(1).evaluate(0);
disc = expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd);
end

function disc = pathcon(expr, N, args)
disc = [];
t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;
for n=1:N
    tt = args.t(n).y(1);
    xx = args.x(n).y(:,1);
    zz = args.z(n).evaluate(0); % Does not have a parameter at tau==0
    uu = args.u(n).y(:);
    dd = args.ders(n).evaluate(0);
    disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
end
tt = args.t(N+1).y;
xx = args.x(N+1).y(:); % Only a point
zz = args.z(N).evaluate(1); % Does not have a parameter at N+1.
uu = args.u(N).y(:); % Same control input as N,
pp = args.p;
tp = args.tps;
int = args.ints;
dd = args.ders(n).evaluate(1);
disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
end

function disc = hard_pathcon(expr, N, args)
disc = [];
t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;
for n=1:N
    tt = args.t(n).y(1);
    xx = args.x(n).y(:,1);
    zz = args.z(n).evaluate(0); % Does not have a parameter at tau==0
    uu = args.u(n).y(:);
    dd = args.ders(n).evaluate(0);
    disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
    for r = 2:length(args.tau)
        tt = args.t(n).y(r);
        xx = args.x(n).y(:, r);
        zz = args.z(n).y(:, r-1);
        uu = args.u(n).y(:);
        dd = args.ders(n).evaluate(args.tau(r));
        disc = [disc, ...
            expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
    end
end
tt = args.t(N+1).y;
xx = args.x(N+1).y(:); % Only a point
zz = args.z(N).evaluate(1); % Does not have a parameter at N+1.
uu = args.u(N).y(:); % Same control input as N,
dd = args.ders(n).evaluate(1);
disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
end

function disc = ival_pathcon(expr, N, args)
disc = [];

I0 = expr.t0;
If = expr.tf;
T0 = args.t(1).t0;
Tf = args.t(1).tf;
I0 = yop.IF(I0 == yop.initial_timepoint, T0, I0);
I0 = yop.IF(I0 == yop.final_timepoint  , Tf, I0);
If = yop.IF(If == yop.initial_timepoint, T0, If);
If = yop.IF(If == yop.final_timepoint  , Tf, If);

t0 = args.t0;
tf = args.tf;
pp = args.p;
tp = args.tps;
int = args.ints;

% Evaluate at the beginning of the interval
e0 = expr.fn( ...
    t0, ...
    tf, ...
    args.t.value(I0), ...
    args.x.value(I0), ...
    args.z.value(I0), ...
    args.u.value(I0), ...
    pp, ...
    tp, ...
    int, ...
    args.ders.value(I0) ...
    );
disc = [disc; e0];

% Evaluate all points within the interval
[n0, r0, nf, rf] = ...
    get_ival_idx(I0, If, args.t(1).t0, args.t(1).tf, N, args.tau);

for n = n0 : min(nf, N) % N+1 is handled by evaluating at If
    for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, length(args.tau))
        tt = args.t(n).y(r);
        xx = args.x(n).y(:, r);
        zz = args.z(n).evaluate(args.tau(r));
        uu = args.u(n).y(:);
        dd = args.ders(n).evaluate(args.tau(r));
        disc = [disc; expr.fn(t0, tf, tt, xx, zz, uu, pp, tp, int, dd)];
    end
end

% Evaluate final point
ef = expr.fn( ...
    t0, ...
    tf, ...
    args.t.value(If), ...
    args.x.value(If), ...
    args.z.value(If), ...
    args.u.value(If), ...
    pp, ...
    tp, ...
    int, ...
    args.ders.value(If));
disc = [disc; ef];
end

function [w_lb, w_ub] = box_bnd(t0, tf, N, tau, ocp)

dt = (tf-t0)/N;
% can be casadi SX/MX, in that case open horizon. If so, it is preferable
% to use a dummy numeric value as it is faster.
dt = yop.IF(isnumeric(dt), dt, 1); 

t0_lb = ocp.t0_lb(t0);
t0_ub = ocp.t0_ub(t0);
tf_lb = ocp.tf_lb(t0);
tf_ub = ocp.tf_ub(t0);

x_ub = ocp.x0_ub(t0);
x_lb = ocp.x0_lb(t0);
t = yop.IF(isnumeric(t0), t0, 1);
for n=1:N
    for r = tau(yop.IF(n==1,2,1):end)        
        x_ub = [x_ub; ocp.x_ub(t+r)];
        x_lb = [x_lb; ocp.x_lb(t+r)];
    end
    t = t + dt;
end
x_ub = [x_ub; ocp.xf_ub(tf)];
x_lb = [x_lb; ocp.xf_lb(tf)];

z_ub = [];
z_lb = [];
t = yop.IF(isnumeric(t0), t0, 1);
for n=1:N
    for r = tau(2:end)
        z_ub = [z_ub; ocp.z_ub(t+r)];
        z_lb = [z_lb; ocp.z_lb(t+r)];
    end
    t = t + dt;
end

u_ub = ocp.u0_ub(t0);
u_lb = ocp.u0_lb(t0);
t = yop.IF(isnumeric(t0), t0, 1);
for n=2:N-1
    u_ub = [u_ub; ocp.u_ub(t)];
    u_lb = [u_lb; ocp.u_lb(t)];
    t = t + dt;
end
u_ub = [u_ub; ocp.uf_ub(tf)];
u_lb = [u_lb; ocp.uf_lb(tf)];

p_lb = ocp.p_lb(t0);
p_ub = ocp.p_ub(t0);

w_ub = vertcat(t0_ub, tf_ub, x_ub, z_ub, u_ub, p_ub);
w_lb = vertcat(t0_lb, tf_lb, x_lb, z_lb, u_lb, p_lb);
end

function [n0, r0, nf, rf] = get_ival_idx(I0, If, T0, Tf, N, tau)
dT = (Tf-T0)/N;

% First point after I0
n0 = 1 + floor((I0-T0)/dT);
if n0 == N+1
    r0 = 1;
else
    t_n0 = dT*(n0-1);
    tau_I0 = (I0-t_n0)/dT;
    % Does not use >= in order to always pick the next point
    r0 = find(tau - tau_I0 > 0, 1); 
    if isempty(r0)
        % Next point is the next control interval
        r0 = 1;
        n0 = n0+1;
    end
end

% Last point before If
nf = 1 + floor((If-T0)/dT);
t_nf = dT*(nf-1);
tau_If = (If-t_nf)/dT;
rf = find(tau - tau_If < 0, 1, 'last');
if isempty(rf)
    if nf==1
        rf = 1;
    else
        rf = length(tau);
        nf = nf-1;
    end
end

end
