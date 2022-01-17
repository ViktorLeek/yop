yops Times: t t0 tf 
yops State: x size: [8,1] scaling: [1e2,1e5,1e5,1e4,1e2,1,10,1e5]
yops Ctrls: u size: [3,1] scaling: [1e2, 1, 1e5] deg: 3
upshift_param;
w_ice=x(1); w_tr=x(5); u_f=u(1); u_wg=u(2); P_gen=u(3);
%% Objective
J = @(dx, y) 1e-4 * int( der(dx(5))^2 ) + 1e-3 * int( y.cylinder.fuel_massflow );

%% First phase
[dx1, y1, c1] = coupled_gear_first(x, u);
p1 = yop.ocp('Up-shift - Phase 1');
p1.min( J(dx1, y1) );
p1.st( t0 == 0 );
p1.st( tf == 0.7 );
p1.st( der(x) == dx1 );
p1.st(  x(t0) == x0 );
p1.st( x_min <= x <= x_max );
p1.st( u_min <= u <= u_max );
p1.st( y1.engine.torque(tf) ==  y1.emachine.torque(tf) );
p1.st( c1{:} );

%% Second phase
[dx2, y2, c2] = decoupled(x, u);
p2 = yop.ocp('Up-shift - Phase 2');
p2.min( J(dx2, y2) );
p2.st( tf-t0 == 0.3 ); % Phase duration
p2.st( der(x) == dx2 );
p2.st(  x_min <=   x   <= x_max  );
p2.st(  u_min <=   u   <= u_max  );
p2.st( y2.engine.torque(tf) == y2.emachine.torque(tf) );
p2.st( w_ice(tf) == w_tr(tf)*i_t(2) ); % Match transmission speed
p2.st( c2{:} );

%% Third phase
[dx3, y3, c3] = coupled_gear_second(x, u);
p3 = yop.ocp('Up-shift - Phase 3');
p3.min( J(dx3, y3) );
p3.st( tf-t0 == 0.504 ); % Phase duration
p3.st( der(x) == dx3 );
p3.st( xf_min <= x(tf) <= xf_max );
p3.st(  x_min <=   x   <= x_max  );
p3.st(  u_min <=   u   <= u_max  );
p3.st( dx3(tf) == 0 );
p3.st( c3{:} );

%% Initial guess - phase 1
e1 = max(0, y1.emachine.torque - y1.engine.torque);
kp1 = 10;
ivp1 = yop.ivp(t0==0, tf==0.7);
ivp1.add( der(x) == dx1 );
ivp1.add(  x(t0) == x0  );
ivp1.add( u_f   == kp1*e1 );
ivp1.add( u_wg  == 1 );
ivp1.add( P_gen == 0 );
sim1 = ivp1.solve();
p1.guess = sim1;

%% Initial guess - phase 2
e2  = w_ice - w_tr*i_t(2);
kp2 = 20e3;
u_pgen = min(kp2*e2, u_max(3));
u_pgen = max(u_pgen, u_min(3));
ivp2 = yop.ivp();
ivp2.add(t0 == sim1.value(tf)      );
ivp2.add(tf == sim1.value(tf) + 0.3);
ivp2.add( der(x) == dx2 );
ivp2.add(  x(t0) == sim1.value(x(tf)) );
ivp2.add( u_f   == 0 );
ivp2.add( u_wg  == 0 );
ivp2.add( P_gen == u_pgen );
sim2 = ivp2.solve();
p2.guess = sim2;

%% Initial guess - phase 3
yops State: I % PI-controller integration state
e3 = xf_max(1) - x(1);
x03 = sim2.value(x(tf));
x03(5) = x03(1)/i_t(2); % should be synchronized
ivp3 = yop.ivp();
ivp3.add(t0 == sim2.value(tf) );
ivp3.add(tf == sim2.value(tf) + 2.0);
ivp3.add( der(x) == dx3 );
ivp3.add(  x(t0) == x03 );
ivp3.add( der(I) == e3 );
ivp3.add(  I(t0) == 0 );
ivp3.add( u_f   == min(max(0, 50*e3 + 20*I), u_max(1)) );
ivp3.add( u_wg  == 0 );
ivp3.add( P_gen == -3*x(8) ); % Returned stored energy
sim3 = ivp3.solve();
p3.guess = sim3;
% sol3 = p3.solve('guess', sim3);

%% Plot initial guess
% figure(1);
% for k=1:8
%     subplot(4,2,k); hold on
%     sim1.plot(t, x(k))
%     sim2.plot(t, x(k))
%     sim3.plot(t, x(k))
% end
% 
% figure(2);
% for k=1:3
%     subplot(3,1,k); hold on
%     sim1.plot(t, u(k))
%     sim2.plot(t, u(k))
%     sim3.plot(t, u(k))
% end

%% OCP
p4 = yop.ocp();
J = p1 + p2 + p3 + p4
% phases = [p1, p2, p3];
% parent = [p1, p1, p2];
% yop.ocp.solve(J, phases, parent);

%% NLP
N = 6;
d = 5;
cp = 'legendre';
nlp1 = yop.direct_collocation(p1.to_canonical(), N, d, cp);
nlp2 = yop.direct_collocation(p2.to_canonical(), N, d, cp);
nlp3 = yop.direct_collocation(p3.to_canonical(), N, d, cp);

%% 1st boundary
% nlps = [nlp1, nlp2];
% nlp2.w_ub(1) =  inf;
% nlp2.w_lb(1) = -inf;
time_continuity = nlp1.tf - nlp2.t0;
state_continuity = nlp1.x(end).eval(0) - nlp2.x(1).eval(0);
param_const = nlp1.p - nlp2.p;
eq = [time_continuity; state_continuity; param_const];
time_positive = nlp2.t0 - nlp2.tf;
g_continuity = [eq; time_positive];
g_continuity_ub = [zeros(size(eq));    0];
g_continuity_lb = [zeros(size(eq)); -inf];

% Merge nlps
nlp = struct;
nlp.J = nlp1.J + nlp2.J;
nlp.w    = [nlp1.w   ; nlp2.w   ];
nlp.w_ub = [nlp1.w_ub; nlp2.w_ub];
nlp.w_lb = [nlp1.w_lb; nlp2.w_lb];
nlp.g    = [nlp1.g   ; nlp2.g   ];
nlp.g_ub = [nlp1.g_ub; nlp2.g_ub];
nlp.g_lb = [nlp1.g_lb; nlp2.g_lb];

% Add continuity constraints
nlp.g    = [nlp.g   ; g_continuity   ];
nlp.g_ub = [nlp.g_ub; g_continuity_ub];
nlp.g_lb = [nlp.g_lb; g_continuity_lb];

%% 2nd boundary
% nlps = [nlps, nlp3];
% nlp3.w_ub(1) =  inf;
% nlp3.w_lb(1) = -inf;
time_continuity = nlp2.tf - nlp3.t0;
state_continuity = nlp2.x(end).eval(0) - nlp3.x(1).eval(0);
param_const = nlp2.p - nlp3.p;
eq = [time_continuity; state_continuity; param_const];
time_positive = nlp3.t0 - nlp3.tf;
g_continuity = [eq; time_positive];
g_continuity_ub = [zeros(size(eq));    0];
g_continuity_lb = [zeros(size(eq)); -inf];

% Merge
nlp.J = nlp.J + nlp3.J;
nlp.w    = [nlp.w   ; nlp3.w];
nlp.w_ub = [nlp.w_ub; nlp3.w_ub];
nlp.w_lb = [nlp.w_lb; nlp3.w_lb];
nlp.g    = [nlp.g   ; nlp3.g];
nlp.g_ub = [nlp.g_ub; nlp3.g_ub];
nlp.g_lb = [nlp.g_lb; nlp3.g_lb];

% Add continuity constraints
nlp.g    = [nlp.g   ; g_continuity   ];
nlp.g_ub = [nlp.g_ub; g_continuity_ub];
nlp.g_lb = [nlp.g_lb; g_continuity_lb];

%% Retriever functions
f1 = casadi.Function('f1', {nlp.w}, {nlp1.w});
f2 = casadi.Function('f2', {nlp.w}, {nlp2.w});
f3 = casadi.Function('f2', {nlp.w}, {nlp3.w});

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

nlp3.ocp_t0 = @(w) nlp3.ocp_t0(full(f3(w)));
nlp3.ocp_tf = @(w) nlp3.ocp_tf(full(f3(w)));
nlp3.ocp_t  = @(w) nlp3.ocp_t (full(f3(w)));
nlp3.ocp_x  = @(w) nlp3.ocp_x (full(f3(w)));
nlp3.ocp_z  = @(w) nlp3.ocp_z (full(f3(w)));
nlp3.ocp_u  = @(w) nlp3.ocp_u (full(f3(w)));
nlp3.ocp_p  = @(w) nlp3.ocp_p (full(f3(w)));

nlps = [nlp1, nlp2, nlp3];
%% NLP Solution
opts.ipopt.max_iter = 10000;
opts.ipopt.acceptable_tol = 1e-7;
solver=casadi.nlpsol('solver','ipopt',struct('f',nlp.J,'x',nlp.w,'g',nlp.g), opts);
nlp_sol = solver( ...
    'x0', [nlp1.w0; nlp2.w0; nlp3.w0], ...
    'lbx', nlp.w_lb, ...
    'ubx', nlp.w_ub, ...
    'ubg', nlp.g_ub, ...
    'lbg', nlp.g_lb ...
    );

%% OCP Solution
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
sol = yop.ocp_sol(mx_vars, unique([p1.ids, p2.ids, p3.ids]), w_opt, 3*N, d, cp);

%% Plot solution
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sol.plot(t, x(k), 'mag', 5)
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sol.plot(t, u(k))
end


%% Initial value
% yops State: I
% x0 = xi_max;
% x0([2,3,6]) = [2e5; 2e5; 0];
% e = 1.466246250000000e+02 - x(1);
% ivp1 = yop.ivp(t0==0, tf==60.5);
% ivp1.add( der(x) == dx1 );
% ivp1.add( der(I) == e );
% ivp1.add(  x(t0) == x0  );
% ivp1.add(  I(t0) == 0 );
% ivp1.add( u_f   == min(max(0, 5*e + 1*I), 150) );
% ivp1.add( u_wg  == 0 );
% ivp1.add( P_gen == 0 );
% sim1 = ivp1.solve();





