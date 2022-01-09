%% Optimal Control Problem Benchmark
yopvar times: t t0 tf
yopvar states: w_ice p_im p_em w_tc scaling: [1e3, 1e5, 1e5, 1e3]
yopvar controls: u_f u_wg P_gen scaling: [1, 1, 1e5]

% States       [rad/s]       [Pa]      [Pa]   [rad/s]
x =     [        w_ice;     p_im;     p_em;     w_tc];
x0    = [rpm2rad( 800); 1.0143e5; 1.0975e5; 2.0502e3]; 
x_min = [rpm2rad( 800);    8.1e4;    9.1e4;      500]; 
x_max = [rpm2rad(2500);    3.5e5;    4.0e5;     15e3]; 

% Control [mg/cyc/cyl]  Effective area [-]      [W]
u =     [        u_f;                u_wg;   P_gen];
u_min = [          0;                   0;       0];
u_max = [        150;                   1;   100e3];

% Dynamics and outputs
[dx, y] = genset_model(x, u);

%% Initial guess
yopvar state: I % PID integral state

% Desired engine speed
wd = rpm2rad(1500);

% PI Engine Speed Controller - with back calculation anti-windup
K  = 10.0; 
Ti =  0.5; 
Tt =  1.0;
e = wd - w_ice;
u_pi = K*e + I;
es = u_f - u_pi; % Non-zero when control is saturated

% Power outtake based on logistic function
P_peak = 100e3; 
T0 = 0.33; 
s = 90;
P_dem = P_peak/(1 + exp(-s*(t-T0)));

sim = yop.ivp(t0==0, tf==1.4);
sim.add( der(x) == dx );
sim.add(  x(t0) == x0 );
sim.add( u_f == smoke_limiter(u_pi, y.u_f_max, u_min(1), u_max(1)) );
sim.add( u_wg == 0 );
sim.add( P_gen == P_dem );
sim.add( der(I) == K/Ti*e + es/Tt ); % PID controller dynamics
sim.add( I(t0)  == 0 );
res = sim.solve('solver', 'ode15s'); 
% res = sim.solve('solver', 'idas', 'points', 100);

%% Plot simulation results
figure(1)
subplot(411); hold on
res.plot(t, rad2rpm(w_ice))
subplot(412); hold on
res.plot(t, p_im)
subplot(413); hold on
res.plot(t, p_em)
subplot(414); hold on
res.plot(t, w_tc)

figure(2)
subplot(311); hold on
res.plot(t, u_f)
res.plot(t, y.u_f_max, '--', 'LineWidth', 2);
subplot(312); hold on
res.stairs(t, u_wg)
subplot(313); hold on
res.plot(t, P_gen)

%% Optimal control problem
ocp = yop.ocp('Optimal Control Problem Benchmark');
ocp.min( 1e3*int(y.cylinder.fuel_massflow) ); % Min fuel mass
% Problem horizon
ocp.st( t0==0, 1.0<=tf<=1.4 );
% Differential constraint
ocp.st( der(x) == dx );
ocp.st(  x(t0) == x0 );
% Box contraints
ocp.st( x_min <= x <= x_max );
ocp.st( u_min <= u <= u_max );
% Path constraints
ocp.st( y.engine.torque >= 0 );
ocp.hard( y.phi <= y.phi_max ); % Constraint applies to all collocation points
% Terminal conditions
ocp.st(  P_gen(tf) == 100e3 ); % [W]
ocp.st( int(P_gen) >= 100e3 ); % [J]
ocp.st(     dx(tf) == 0 );     % Stationarity
    
sol = ocp.solve('intervals', 75, 'degree', 3,'guess', res);

%% Plot the optimal control and trajectory
figure(1)
subplot(411); hold on
sol.plot(t, rad2rpm(w_ice), 'mag', 2)
subplot(412); hold on
sol.plot(t, p_im, 'mag', 2)
subplot(413); hold on
sol.plot(t, p_em, 'mag', 2)
subplot(414); hold on
sol.plot(t, w_tc, 'mag', 2)

figure(2)
subplot(311); hold on
sol.stairs(t, u_f)
sol.plot(t, y.u_f_max);
subplot(312); hold on
sol.stairs(t, u_wg)
subplot(313); hold on
sol.plot(t, P_gen)
