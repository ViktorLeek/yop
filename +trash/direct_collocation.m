function nlp = direct_collocation(ocp, N, d, cp)

% If horizon is not fixed, T0 and Tf are nonsense. But since they are only
% to be used for a fixed horizon, that is OK.
[~, T0, Tf] = fixed_horizon(ocp);

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
    x_n = yop.cx(['x_' num2str(n)], n_x(ocp), d+1);
    x(n) = yop.interpolating_poly(tau, x_n, T0, Tf, N);
end
x_n = yop.cx(['x_' num2str(N+1)], n_x(ocp));
x(N+1) = yop.interpolating_poly(0, x_n, T0, Tf, N);

z = yop.interpolating_poly.empty(N, 0);
for n=1:N
    z_n = yop.cx(['z_' num2str(n)], n_z(ocp), d);
    z(n) = yop.interpolating_poly(tau(2:end), z_n, T0, Tf, N);
end

u = yop.interpolating_poly.empty(N, 0);
for n=1:N
    u_n = yop.cx(['u_' num2str(n)], n_u(ocp));
    u(n) = yop.interpolating_poly(0, u_n, T0, Tf, N);
end

p = yop.cx('p', n_p(ocp));

[tps, ints, ders] = param_special_nodes(ocp,N,tau,dt,t0,tf,t,x,z,u,p,T0,Tf);

J = ocp.objective.fn(t0, tf, p, tps, ints);

g = discretize_ode(ocp,N,tau,dt,t0,tf,t,x,z,u,p,tps,ints,ders);

for pc = ocp.equality_constraints
    disc = parameterize_expr(pc,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders);
    g = [g; disc(:)];
end

h = [];
for pc = ocp.inequality_constraints
    disc = parameterize_expr(pc,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders);
    h = [h; disc(:)];
end

[w_lb, w_ub] = box_bnd(T0, Tf, N, tau, ocp);
w = vertcat(t0, tf, vec(x), vec(z), vec(u), vec(p));
t0 = casadi.Function('t0', {w}, {t0});
tf = casadi.Function('tf', {w}, {tf});
t  = casadi.Function('t' , {w}, {mat(t)});
x  = casadi.Function('x' , {w}, {mat(x)});
z  = casadi.Function('z' , {w}, {mat(z)});
u  = casadi.Function('u' , {w}, {mat(u)});
p  = casadi.Function('p' , {w}, {p});

nlp = struct;
nlp.J = J;
nlp.w = w;
nlp.w_ub = w_ub;
nlp.w_lb = w_lb;
nlp.g = g;
nlp.g_ub = zeros(size(g));
nlp.g_lb = zeros(size(g));
nlp.h = h;
nlp.h_ub = zeros(size(h));
nlp.h_lb = -inf(size(h));
nlp.t0 = @(w) full(t0(w));
nlp.tf = @(w) full(tf(w));
nlp.t  = @(w) full(t(w));
nlp.x  = @(w) full(x(w));
nlp.z  = @(w) full(z(w));
nlp.u  = @(w) full(u(w));
nlp.p  = @(w) full(p(w));
end


function [tps, ints, ders] = param_special_nodes ...
    (ocp, N, tau, dt, t0, tf, t, x, z, u, p, T0, Tf)

tps = [];
ints = [];
ders = yop.interpolating_poly.empty(N, 0);
for n=1:N
    ders(n) = yop.interpolating_poly(tau, [], T0, Tf, N);
end

for node = ocp.special_nodes
    tmp_tp  = [tps;  zeros(ocp.n_tp  - length(tps), 1)];
    tmp_int = [ints; zeros(ocp.n_int - length(ints), 1)];
    switch node.type
        case yop.ocp_expr.tp_type
            tp = parameterize_timepoint(node, t0, tf, t, x, z, ...
                u, p, tmp_tp, tmp_int, ders, ocp.n_der);
            tps = [tps; tp(:)];
            
        case yop.ocp_expr.int_type
            int = parameterize_integral(node, N, tau, dt, t0, tf, t, x, ...
                z, u, p, tmp_tp, tmp_int, ders, ocp.n_der);
            ints = [ints; int(:)];
            
        case yop.ocp_expr.der_type
            ders = parameterize_derivative(node, N, tau, dt, t0, tf, t, ...
                x, z, u, p, tmp_tp, tmp_int, ders, ocp.n_der);
            
        otherwise
            error(yop.msg.unexpected_error);
    end
end
end

function val = parameterize_timepoint(tp,t0,tf,t,x,z,u,p,tps,ints,ders,n_der)
tt = t.value(tp.timepoint);
xx = x.value(tp.timepoint);
zz = z.value(tp.timepoint);
uu = u.value(tp.timepoint); 
dd = ders.value(tp.timepoint);
dd = [dd;  zeros(n_der  - length(dd), 1)];
val = tp.fn(t0, tf, tt, xx, zz, uu, p, tps, ints, dd);
end

function I=parameterize_integral(i,N,tau,dt,t0,tf,t,x,z,u,p,tps,ints,ders,n_der)
I = 0;
for n=1:N
    yval = [];
    for r=1:length(tau)
        tt = t(n).y(r);
        xx = x(n).y(:, r);
        zz = z(n).evaluate(tau(r));
        uu = u(n).y;
        pp = p;
        dd = ders(n).evaluate(tau(r));
        dd = [dd;  zeros(n_der  - length(dd), 1)];
        val_r = i.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd);
        yval = [yval, val_r(:)];
    end
    lp = yop.lagrange_polynomial(tau, yval).integrate();
    I = I + lp.evaluate(1)*dt;
end
end

function ders = parameterize_derivative ...
    (der, N, tau, dt, t0, tf, t, x, z, u, p, tps, ints, ders, n_der)

for n=1:N
    yval = [];
    for r=1:length(tau)
        tt = t(n).y(r);
        xx = x(n).y(:, r);
        zz = z(n).evaluate(tau(r));
        uu = u(n).y;
        pp = p;
        dd = ders(n).evaluate(tau(r));
        dd = [dd;  zeros(n_der  - length(dd), 1)];
        val_r = der.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd);
        yval = [yval, val_r(:)];
    end
    yn = yop.lagrange_polynomial(tau, ...
        yval).differentiate().evaluate(tau)/dt;
    ders(n).y = [ders(n).y; yn];
end
end

function g = discretize_ode(ocp,N,tau,dt,t0,tf,t,x,z,u,p,tps,ints,ders)
g = []; % Equality constraints from discretization
for n=1:N % Dynamics
    dx = x(n).differentiate();
    for r=2:length(tau)
        dxr = dx.evaluate(tau(r));
        tt = t(n).y(r);
        xx = x(n).y(:, r);
        zz = z(n).y(:, r-1); % only has d parameters
        uu = u(n).y(:);
        pp = p;
        dd = ders(n).evaluate(tau(r));
        f = ocp.ode.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd);
        a = ocp.alg.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd);
        g = [g; (dxr - dt*f); a];
    end
end

for n=1:N % Continuity
    gn = x(n).evaluate(1) - x(n+1).evaluate(0);
    g = [g; gn];
end
end

function disc = parameterize_expr(expr,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders)
if is_transcription_invariant(expr)
    disc = parameterize_invariant(expr,t0,tf,t,x,z,u,p,tps,ints,ders);
elseif is_ival(expr)
    disc = parameterize_ival(expr,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders);
else
    disc = parameterize_all(expr,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders);
end
end

function disc = parameterize_invariant(expr,t0,tf,t,x,z,u,p,tps,ints,ders)
tt = t(1).evaluate(0);
xx = x(1).evaluate(0);
zz = z(1).evaluate(0);
uu = u(1).evaluate(0);
dd = ders(1).evaluate(0);
disc = expr.fn(t0, tf, tt, xx, zz, uu, p, tps, ints, dd);
end

function disc = parameterize_all(expr,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders)
disc = [];
for n=1:N
    tt = t(n).y(1);
    xx = x(n).y(:,1);
    zz = z(n).evaluate(0); % Does not have a parameter at tau==0
    uu = u(n).y(:);
    pp = p;
    dd = ders(n).evaluate(0);
    disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd)];
    if expr.is_hard
        for r = 2:length(tau)
            tt = t(n).y(r);
            xx = x(n).y(:, r);
            zz = z(n).y(:, r-1);
            uu = u(n).y(:);
            dd = ders(n).evaluate(tau(r));
            disc = [disc, ...
                expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd)];
        end
    end
end
tt = t(N+1).y;
xx = x(N+1).y(:); % Only a point
zz = z(N).evaluate(1); % Does not have a parameter at N+1.
uu = u(N).y(:); % Same control input as N,
pp = p;
dd = ders(n).evaluate(1);
disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, dd)];
end

function disc = parameterize_ival(expr,N,tau,t0,tf,t,x,z,u,p,tps,ints,ders)
disc = [];

[I0, If] = get_ival(expr);
T0 = t(1).t0;
Tf = t(1).tf;
I0 = yop.IF(I0==yop.initial_timepoint, T0, I0);
I0 = yop.IF(I0==yop.final_timepoint  , Tf, I0);
If = yop.IF(If==yop.initial_timepoint, T0, If);
If = yop.IF(If==yop.final_timepoint  , Tf, If);

% Evaluate at the beginning of the interval
disc = [disc, expr.fn(t0, tf, t.value(I0), x.value(I0), z.value(I0), ...
    u.value(I0), p, tps, ints, ders.value(I0))];

% Evaluate all points within the interval
[n0, r0, nf, rf] = get_ival_idx(I0, If, t(1).t0, t(1).tf, N, tau);
for n=n0:min(nf, N) % N+1 is handled by evaluating at If
    % Intervals are treated as hard constraints.
    for r = yop.IF(n==n0, r0, 1) : yop.IF(n==nf, rf, length(tau))
        tt = t(n).y(r);
        xx = x(n).y(:, r);
        zz = z(n).evaluate(tau(r));
        uu = u(n).y(:);
        dd = ders(n).evaluate(tau(r));
        disc = [disc, ...
            expr.fn(t0, tf, tt, xx, zz, uu, p, tps, ints, dd)];
    end
end

% Evaluate final point
disc = [disc, expr.fn(t0, tf, t.value(If), x.value(If), z.value(If), ...
    u.value(If), p, tps, ints, ders.value(If))];
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

function [w_lb, w_ub] = box_bnd(t0, tf, N, tau, ocp)

dt = (tf-t0)/N;
% can be casadi SX/MX, in that case open horizon. If so, it is preferable
% to use a dummy numeric value as it is a lot faster.
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

%%

% t0_lb = ocp.t0_lb;
% t0_ub = ocp.t0_ub;
% tf_lb = ocp.tf_lb;
% tf_ub = ocp.tf_ub;
% 
% reps = N*length(tau) + 1;
% x_lb = repmat(ocp.x_lb(1), reps, 1);
% x_ub = repmat(ocp.x_ub(1), reps, 1);
% x_lb(1 : ocp.n_x) = ocp.x0_lb(1);
% x_ub(1 : ocp.n_x) = ocp.x0_ub(1);
% x_lb(end - ocp.n_x + 1 : end) = ocp.xf_lb(1);
% x_ub(end - ocp.n_x + 1 : end) = ocp.xf_ub(1);
% 
% z_lb = repmat(ocp.z_lb(1), N*(length(tau)-1), 1);
% z_ub = repmat(ocp.z_ub(1), N*(length(tau)-1), 1);
% 
% u_lb = repmat(ocp.u_lb(1), N, 1);
% u_ub = repmat(ocp.u_ub(1), N, 1);
% u_lb(1 : ocp.n_u) = ocp.u0_lb(1);
% u_ub(1 : ocp.n_u) = ocp.u0_ub(1);
% u_lb(end - ocp.n_u + 1 : end) = ocp.uf_lb(1);
% u_ub(end - ocp.n_u + 1 : end) = ocp.uf_ub(1);
% 
% p_lb = ocp.p_lb(1);
% p_ub = ocp.p_ub(1);

w_ub = vertcat(t0_ub, tf_ub, x_ub, z_ub, u_ub, p_ub);
w_lb = vertcat(t0_lb, tf_lb, x_lb, z_lb, u_lb, p_lb);
end
