%% Multp-phase Bryson-Denham Problem
yops times: t t0 tf % Parsed by position: t, t0, tf
yops states: x v    % position, speed
yops ctrls: a       % acceleration

p1 = yop.ocp('Bryson - Phase 1');
p1.min( 1/2 * int(a^2) );
p1.st( t0 == 0 );
p1.st( tf == 0.5 );
p1.st( der(x) == v );
p1.st( der(v) == a );
p1.st( x(t0)  == 0 );
p1.st( v(t0)  == 1 );
p1.st( x <= 1/9 );

p2 = yop.ocp('Bryson - Phase 2');
p2.min( 1/2 * int(a^2) );
p2.st( tf == 1 );
p2.st( der(x) == v );
p2.st( der(v) == a );
p2.st( x(tf)  == 0 );
p2.st( v(tf)  == -1 );
p2.st( x <= 1/9 );

%% NLP
N = 10;
d = 2;
cp = 'legendre';
nlp1 = yop.direct_collocation(p1.to_canonical(), N, d, cp);
nlp2 = yop.direct_collocation(p2.to_canonical(), N, d, cp);

% Continuity constratints
% nlp2.w_ub(1) =  inf;
% nlp2.w_lb(1) = -inf;
time_continuity  = nlp1.tf - nlp2.t0;
time_positive    = nlp2.t0 - nlp2.tf;
state_continuity = nlp1.x(end).eval(0) - nlp2.x(1).eval(0);
parameter_const  = nlp1.p - nlp2.p;
eq               = [time_continuity; state_continuity; parameter_const];
g_continuity     = [eq; time_positive];
g_continuity_ub  = [zeros(size(eq));    0];
g_continuity_lb  = [zeros(size(eq)); -inf];

% Merge nlps
nlp = struct;
nlp.J = nlp1.J + nlp2.J;
nlp.w    = [nlp1.w   ; nlp2.w   ];
nlp.w_ub = [nlp1.w_ub; nlp2.w_ub];
nlp.w_lb = [nlp1.w_lb; nlp2.w_lb];
nlp.g    = [nlp1.g   ; nlp2.g   ];
nlp.g_ub = [nlp1.g_ub; nlp2.g_ub];
nlp.g_lb = [nlp1.g_lb; nlp2.g_lb];

% Add continuoty constraints
nlp.g    = [nlp.g   ; g_continuity   ];
nlp.g_ub = [nlp.g_ub; g_continuity_ub];
nlp.g_lb = [nlp.g_lb; g_continuity_lb];

f1 = casadi.Function('f1', {nlp.w}, {nlp1.w});
f2 = casadi.Function('f2', {nlp.w}, {nlp2.w});

nlp1.ocp_t0 = @(w) nlp1.ocp_t0(full(f1(w)));
nlp1.ocp_tf = @(w) nlp1.ocp_tf(full(f1(w)));
nlp1.ocp_t  = @(w) nlp1.ocp_t (full(f1(w)));
nlp1.ocp_x  = @(w) nlp1.ocp_x (full(f1(w)));
nlp1.ocp_z  = @(w) nlp1.ocp_z (full(f1(w)));
nlp1.ocp_u  = @(w) nlp1.ocp_u (full(f1(w)));
nlp1.ocp_p  = @(w) nlp1.ocp_p (full(f1(w)));

nlp2.ocp_t0 = @(w) nlp2.ocp_t0(full(f2(w)));
nlp2.ocp_tf = @(w) nlp2.ocp_tf(full(f2(w)));
nlp2.ocp_t  = @(w) nlp2.ocp_t (full(f2(w)));
nlp2.ocp_x  = @(w) nlp2.ocp_x (full(f2(w)));
nlp2.ocp_z  = @(w) nlp2.ocp_z (full(f2(w)));
nlp2.ocp_u  = @(w) nlp2.ocp_u (full(f2(w)));
nlp2.ocp_p  = @(w) nlp2.ocp_p (full(f2(w)));

solver = casadi.nlpsol('solver', 'ipopt', ...
    struct('f', nlp.J, 'x', nlp.w, 'g', nlp.g));

nlp_sol = solver( ...
    'x0', ones(size(nlp.w)), ...
    'lbx', nlp.w_lb, ...
    'ubx', nlp.w_ub, ...
    'ubg', nlp.g_ub, ...
    'lbg', nlp.g_lb ...
    );

nlps = [nlp1, nlp2];
tt0 = nlps(1).ocp_t0(nlp_sol.x);
tt = nlps(1).ocp_t(nlp_sol.x);
xx = nlps(1).ocp_x(nlp_sol.x);
zz = nlps(1).ocp_z(nlp_sol.x);
uu = nlps(1).ocp_u(nlp_sol.x);
pp = nlps(1).ocp_p(nlp_sol.x);
for nlp_k = nlps(2:end)
    t_k = nlp_k.ocp_t(nlp_sol.x);
    x_k = nlp_k.ocp_x(nlp_sol.x);
    z_k = nlp_k.ocp_z(nlp_sol.x);
    u_k = nlp_k.ocp_u(nlp_sol.x);
    tt = [tt, t_k(2:end)];
    xx = [xx, x_k(:, 2:end)];
    zz = [zz, z_k];
    uu = [uu, u_k];
end
ttf = nlps(end).ocp_tf(nlp_sol.x);

w_opt = struct;
w_opt.t0 = p1.descale_t0(tt0);
w_opt.tf = p1.descale_tf(ttf);
w_opt.t  = p1.descale_t (tt);
w_opt.x  = p1.descale_x (xx);
w_opt.z  = p1.descale_z (zz);
w_opt.u  = p1.descale_u (uu);
w_opt.p  = p1.descale_p (pp);
mx_vars = p1.mx_vars();
sol = yop.ocp_sol(mx_vars, unique([p1.ids, p2.ids]), w_opt, 2*N, d, cp);

w_opt1 = struct;
w_opt1.t0 = p1.descale_t0(nlp1.ocp_t0(nlp_sol.x));
w_opt1.tf = p1.descale_tf(nlp1.ocp_tf(nlp_sol.x));
w_opt1.t  = p1.descale_t (nlp1.ocp_t (nlp_sol.x));
w_opt1.x  = p1.descale_x (nlp1.ocp_x (nlp_sol.x));
w_opt1.z  = p1.descale_z (nlp1.ocp_z (nlp_sol.x));
w_opt1.u  = p1.descale_u (nlp1.ocp_u (nlp_sol.x));
w_opt1.p  = p1.descale_p (nlp1.ocp_p (nlp_sol.x));
sol1 = yop.ocp_sol(mx_vars, p1.ids, w_opt1, N, d, cp);

w_opt2 = struct;
w_opt2.t0 = p2.descale_t0(nlp2.ocp_t0(nlp_sol.x));
w_opt2.tf = p2.descale_tf(nlp2.ocp_tf(nlp_sol.x));
w_opt2.t  = p2.descale_t (nlp2.ocp_t (nlp_sol.x));
w_opt2.x  = p2.descale_x (nlp2.ocp_x (nlp_sol.x));
w_opt2.z  = p2.descale_z (nlp2.ocp_z (nlp_sol.x));
w_opt2.u  = p2.descale_u (nlp2.ocp_u (nlp_sol.x));
w_opt2.p  = p2.descale_p (nlp2.ocp_p (nlp_sol.x));
sol2 = yop.ocp_sol(mx_vars, p2.ids, w_opt2, N, d, cp);

figure(1);
subplot(311); hold on
sol.plot(t, x);
sol1.plot(t, x);
sol2.plot(t0, x(t0), 'o');
subplot(312); hold on
sol.plot(t, v);
sol1.plot(t, v);
sol2.plot(t0, v(t0), 'o');
subplot(313); hold on
sol.stairs(t, a);
sol1.stairs(t, a);
sol2.stairs(t0, a(t0), 'o');
%%
w_opt = struct;
w_opt.t0 = ocp.descale_t0(nlp.t0(nlp_sol.x));
w_opt.tf = ocp.descale_tf(nlp.tf(nlp_sol.x));
w_opt.t  = ocp.descale_t (nlp.t (nlp_sol.x));
w_opt.x  = ocp.descale_x (nlp.x (nlp_sol.x));
w_opt.z  = ocp.descale_z (nlp.z (nlp_sol.x));
w_opt.u  = ocp.descale_u (nlp.u (nlp_sol.x));
w_opt.p  = ocp.descale_p (nlp.p (nlp_sol.x));

sol = yop.ocp_sol(ocp.mx_vars(), ocp.ids, w_opt, N, d, cp);
%%

ocp = p1 + p2; % Sum the objectives
[sol, sol1, sol2] = ocp.solve('intervals', [25, 25]);

figure(1);
subplot(311); hold on
sol.plot(t, x);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.stairs(t, a);