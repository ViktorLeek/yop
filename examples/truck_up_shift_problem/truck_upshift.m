yops Times: t t0 tf 
yops State: x size: [8,1] scaling: [1e2,1e5,1e5,1e4,1e2,1,10,1e5]
yops Ctrls: u_f   scaling: 1e2         % pw constant
yops Ctrls: u_wg  scaling:   1 deg: 0  % pw linear
yops Ctrls: P_gen scaling: 1e5 deg: 0  % pw linear
% P_gen.du.scaling(1e7); % Scale derivative of P_gen
u = [u_f; u_wg; P_gen];
w_ice=x(1); w_tr=x(5);
upshift_param;

%% First phase
[dx1, y1, c1] = coupled_gear_first(x, u);
% J1 = 1e-5 * 1/(tf-t0)*int( der( dx1(5)^2 ) ); % Minimize driveline jerk;
% J1 = 1e-5 * int( der( dx1(5)^2 ) ); % Minimize driveline jerk;
J1 = int( y1.cylinder.fuel_massflow ); % Minimize driveline jerk;
p1 = yop.ocp('Up-shift - Phase 1');
p1.min( J1 );
p1.st( t0 == 0 );
p1.st( tf == 0.7 );
p1.st( der(x)  == dx1 );
% p1.st( dx1(t0) == 0 );
p1.st( P_gen(t0) == 0 );
p1.st( xi_min <= x(t0) <= xi_max ); 
p1.st(  x_min <=  x(t) <= x_max  );
p1.st(  u_min <=   u   <= u_max  );
% p1.st(   -3 <= der(u_wg) <= 3    );
% p1.st( -1e4 <= der(y1.generator.torque) <= 1e4 );
p1.st( y1.engine.torque(tf) ==  y1.emachine.torque(tf) );
p1.st( c1{:} );

%% Second phase
[dx2, y2, c2] = decoupled(x, u);
% J2 = 1e-8 * int( der( dx2(5)^2 ) ); % Minimize driveline jerk
J2 = int( y2.cylinder.fuel_massflow ); % Minimize driveline jerk;
p2 = yop.ocp('Up-shift - Phase 2');
p2.min( J2 );
p2.st( tf-t0 == 0.3 ); % Phase duration
p2.st( der(x) == dx2 );
p2.st(  x_min <=   x   <= x_max  );
p2.st(  u_min <=   u   <= u_max  );
% p2.st(   -3 <= der(u_wg) <= 3    );
% p2.st( -1e4 <= der(y2.generator.torque) <= 1e4 );
p2.st( y2.engine.torque(tf) == y2.emachine.torque(tf) );
p2.st( w_ice(tf) == w_tr(tf)*i_t(2) ); % Match transmission speed
p2.st( c2{:} );

%% Third phase
[dx3, y3, c3] = coupled_gear_second(x, u);
% J3 = 1e-8 * int( der( dx3(5)^2 ) ); % Minimize driveline jerk
J3 = int( y3.cylinder.fuel_massflow ); % Minimize driveline jerk;
p3 = yop.ocp('Up-shift - Phase 3');
p3.min( J3 );
p3.st( tf-t0 == 0.504 ); % Phase duration
p3.st( der(x) == dx3 );
p3.st( xf_min <= x(tf) <= xf_max );
p3.st(  x_min <=   x   <= x_max  );
p3.st(  u_min <=   u   <= u_max  );
% p3.st(   -3 <= der(u_wg) <= 3    );
% p3.st( -1e4 <= der(y3.generator.torque) <= 1e4 );
% p3.st( dx3(tf) == 0 );
p3.st( c3{:} );

%% Initial guess - phase 1 
% ------------------------------------- Fixar ett stationärt
% begynnelsevärde
yops State: I
x0 = xi_max;
x0([2,3,6]) = [2e5; 2e5; 0];
e1 = xi_min(7) - x(7);
kp1 = 10;
ivp1 = yop.ivp(t0==0, tf==10.5);
ivp1.add( der(x) == dx1 );
ivp1.add(  x(t0) == x0  );
ivp3.add( der(I) == e1 );
ivp3.add(  I(t0) == 0 );
ivp1.add( u_f   == kp1*e1 );
ivp1.add( u_wg  == 0 );
ivp1.add( P_gen == 0 );
sim1 = ivp1.solve();
p1.guess = sim1;

%% Initial guess - phase 1
x0 = xi_max;
x0([2,3,6]) = [2e5; 2e5; 0];
e1 = max(0, y1.emachine.torque - y1.engine.torque);
kp1 = 10;
ivp1 = yop.ivp(t0==0, tf==10.5);
ivp1.add( der(x) == dx1 );
ivp1.add(  x(t0) == x0  );
ivp1.add( u_f   == kp1*e1 );
ivp1.add( u_wg  == 0 );
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

%% Plot initial guess
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sim1.plot(t, x(k))
%     sim2.plot(t, x(k))
%     sim3.plot(t, x(k))
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sim1.plot(t, u(k))
%     sim2.plot(t, u(k))
%     sim3.plot(t, u(k))
end

%% NLP
N = 50;
d = 4;
cp = 'legendre';
nlp1 = yop.direct_collocation(p1.to_canonical(), N, d, cp);
nlp2 = yop.direct_collocation(p2.to_canonical(), N, d, cp);
nlp3 = yop.direct_collocation(p3.to_canonical(), N, d, cp);

%% 1st boundary
nlps = [nlp1, nlp2];
nlp2.w_ub(1) =  inf;
nlp2.w_lb(1) = -inf;
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
nlps = [nlps, nlp3];
nlp3.w_ub(1) =  inf;
nlp3.w_lb(1) = -inf;
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

%% Solution
solver=casadi.nlpsol('solver','ipopt',struct('f',nlp.J,'x',nlp.w,'g',nlp.g));
nlp_sol = solver( ...
    'x0', [nlp1.w0; nlp2.w0; nlp3.w0], ...
    'lbx', nlp.w_lb, ...
    'ubx', nlp.w_ub, ...
    'ubg', nlp.g_ub, ...
    'lbg', nlp.g_lb ...
    );












