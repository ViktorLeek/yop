%% Optimal Control Benchmark - with bootstrapping
% Sometimes it is difficult to obtain a good initial guess. For such cases
% a "bootstrapping" strategy can be a good option. In such a scenario the
% problem is relaxed and solved using without an initial guess or a very
% crude one. That solution can then be fed into a less relaxed formulation
% until eventually a guess for the full problem is obtained, and the
% problem can be solved. 
% 
% Here this principle is demonstrated by first using a static guess and a
% relaxed problem (stationarity is the final bound is omitted) and then
% that solution is used to initialize the full problem. This problem is a
% little too simple for this method to be really useful as it is solvable
% without it, but it demonstrates the method.
%% Variables, bounds, dynamics 
yops Times: t t0 tf
yops States: w_ice p_im p_em w_tc nominal: [1e3, 1e5, 1e5, 1e3]
yops Controls: u_f u_wg P_gen     nominal: [1, 1, 1e5]

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
% Notice that x is a [4,1] vector. yop.guess assumes a [length(t_guess),4] 
% vector. Here t is not specified, so yop expectes the guess for x to have 
% size [1, 4] (constant) or [2, 4] (boundary values).
xg = [rpm2rad(1500); 1.3e5; 1.3e5; rpm2rad(50e3)];
ug = [50; 0; 1e3];
guess = yop.guess(t0, 0, tf, 1, x, xg, u, ug);

%% Optimal control problem - Relaxed problem
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
ocp.hard( y.phi <= y.phi_max );
ocp.st( int(P_gen) == 100e3 ); % [J]
% Terminal conditions
ocp.st(  P_gen(tf) == 100e3 ); % [W]
sol = ocp.solve('ival', 25, 'dx', 2,'guess', guess);

%% Full problem - include stationarity constraint and solution refinement
ocp.st( dx(tf) == 0 ); % Stationarity
sol = ocp.solve('ival', 75, 'dx', 3,'guess', sol);

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
