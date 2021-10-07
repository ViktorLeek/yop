function nlp = direct_collocation(ocp, N, d, cp)

% If horizon is not fixed, T0 and Tf are nonsense. But since they are only
% to be used for a fixed horizon, that is OK.
[~, T0, Tf] = fixed_horizon(ocp);

t0 = yop.cx('t0');
tf = yop.cx('tf');
dt = (tf - t0)/N;
tau = [0, casadi.collocation_points(d, cp)];

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

[tps, ints, ders] = param_special_nodes(ocp,N,tau,dt,t0,tf,t,x,z,u,p);

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

w = vertcat(t0, tf, vec(x), vec(z), vec(u), vec(p));
[w_lb, w_ub] = box_bnd(N, d, ocp);

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
nlp.t0 = casadi.Function('t0', {w}, {t0});
nlp.tf = casadi.Function('tf', {w}, {tf});
nlp.t  = casadi.Function('t' , {w}, {mat(t)});
nlp.x  = casadi.Function('x' , {w}, {mat(x)});
nlp.z  = casadi.Function('z' , {w}, {mat(z)});
nlp.u  = casadi.Function('u' , {w}, {mat(u)});
nlp.p  = casadi.Function('p' , {w}, {p});
end


function [tps, ints, ders] = param_special_nodes ...
    (ocp, N, tau, dt, t0, tf, t, x, z, u, p)

tps = [];
ints = [];
ders = [];
for node = ocp.special_nodes
    tmp_tp  = [tps;  zeros(ocp.n_tp  - length(tps), 1)];
    tmp_int = [ints; zeros(ocp.n_int - length(ints), 1)];
    switch node.type
        case yop.ocp_expr.tp
            tp = parameterize_timepoint(node, t0, tf, t, x, z, ...
                u, p, tmp_tp, tmp_int, ders);
            tps = [tps; tp(:)];
            
        case yop.ocp_expr.int
            int = parameterize_integral(node, N, tau, dt, t0, tf, t, x, ...
                z, u, p, tmp_tp, tmp_int, ders);
            ints = [ints; int(:)];
            
        case yop.ocp_expr.der
            error(yop.msg.not_implemented);
            warning(['When implementing dont forget to ' ...
                'multiply derivative with step length']);
            
        otherwise
            error(yop.msg.unexpected_error);
    end
end
end

function val = parameterize_timepoint(tp,t0,tf,t,x,z,u,p,tps,ints,ders)
tt = t.value(tp.timepoint);
xx = x.value(tp.timepoint);
zz = z.value(tp.timepoint);
uu = u.value(tp.timepoint); 
val = tp.fn(t0, tf, tt, xx, zz, uu, p, tps, ints, ders);
end

function I=parameterize_integral(i,N,tau,dt,t0,tf,t,x,z,u,p,tps,ints,ders)

I = 0;
for n=1:N
    yval = [];
    for r=1:length(tau)
        tt = t(n).y(r);
        xx = x(n).y(:, r);
        zz = z(n).evaluate(tau(r)); % Does not have a parameter at tau==0
        uu = u(n).y(:);
        pp = p;
        val_r = i.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders);
        yval = [yval, val_r(:)];
    end
    lp = yop.lagrange_polynomial(tau, yval).integrate();
    I = I + lp.evaluate(1)*dt;
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
        f = ocp.ode.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders);
        a = ocp.alg.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders);
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
    tt = t(1).evaluate(0);
    xx = x(1).evaluate(0);
    zz = z(1).evaluate(0);
    uu = u(1).evaluate(0);
    % dd = der(1).evaluate(0);
    disc = expr.fn(t0, tf, tt, xx, zz, uu, p, tps, ints, ders);
else
    disc = [];
    for n=1:N
        tt = t(n).y(1);
        xx = x(n).y(:,1);
        zz = z(n).evaluate(0); % Does not have a parameter at tau==0
        uu = u(n).y(:);
        pp = p;
        % dd = der(n).evaluate(0);
        disc = [disc, ...
            expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders)];
        if expr.is_hard
            for r = 2:length(tau)
                tt = t(n).y(r);
                xx = x(n).y(:, r);
                zz = z(n).y(:, r-1);
                uu = u(n).y(:);
                disc = [disc, ...
                    expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders)];
            end
        end
    end
    tt = t(N+1).y(1);
    xx = x(N+1).y(:); % Only a point
    zz = z(N).evaluate(1); % Does not have a parameter at N+1.
    uu = u(N).y(:); % Same control input as N, 
    pp = p;
    disc = [disc, expr.fn(t0, tf, tt, xx, zz, uu, pp, tps, ints, ders)];
end
end

function [w_lb, w_ub] = box_bnd(N, d, ocp)
t0_lb = ocp.t0_lb;
t0_ub = ocp.t0_ub;
tf_lb = ocp.tf_lb;
tf_ub = ocp.tf_ub;

reps = N*(d + 1) + 1;
x_lb = repmat(ocp.x_lb, reps, 1);
x_ub = repmat(ocp.x_ub, reps, 1);
x_lb(1 : ocp.n_x) = ocp.x0_lb;
x_ub(1 : ocp.n_x) = ocp.x0_ub;
x_lb(end - ocp.n_x + 1 : end) = ocp.xf_lb;
x_ub(end - ocp.n_x + 1 : end) = ocp.xf_ub;

z_lb = repmat(ocp.z_lb, N*d, 1);
z_ub = repmat(ocp.z_ub, N*d, 1);

u_lb = repmat(ocp.u_lb, N, 1);
u_ub = repmat(ocp.u_ub, N, 1);
u_lb(1 : ocp.n_u) = ocp.u0_lb;
u_ub(1 : ocp.n_u) = ocp.u0_ub;
u_lb(end - ocp.n_u + 1 : end) = ocp.uf_lb;
u_ub(end - ocp.n_u + 1 : end) = ocp.uf_ub;

p_lb = ocp.p_lb;
p_ub = ocp.p_ub;

w_ub = vertcat(t0_ub, tf_ub, x_ub, z_lb, u_ub, p_ub);
w_lb = vertcat(t0_lb, tf_lb, x_lb, z_ub, u_lb, p_lb);
end
