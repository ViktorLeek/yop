%% Optimal Control Problem Benchmark
t0    = yop.time0('t0');
tf    = yop.timef('tf');
t     = yop.time('t');

w_ice = yop.state('w_ice'); % Engine angular velocity
p_im  = yop.state('p_im');  % Intake manifold pressure
p_em  = yop.state('p_em');  % Exhause manifold pressure
w_tc  = yop.state('w_tc');  % Turbocharger angular velocity

u_f = yop.control('name', 'u_f', 'pw', 'linear');   % Fuel injection per cycle per cylinder
u_wg  = yop.control('name', 'u_wg');  % Wastegate control 0-close, 1-fully open
P_gen = yop.control('name', 'P_gen', 'pw', 'linear'); % Generator power

%             [rad/s]       [Pa]      [Pa]   [rad/s]
x =     [        w_ice;     p_im;     p_em;     w_tc]; % State vector
x0    = [rpm2rad( 800); 1.0143e5; 1.0975e5; 2.0502e3];
x_min = [rpm2rad( 800);    8.1e4;    9.1e4;      500];
x_max = [rpm2rad(2500);    3.5e5;    4.0e5;     15e3];

%        [mg/cyc/cyl]  Effective area [-]      [W]
u =     [        u_f;                u_wg;   P_gen]; % Control vector
u_min = [          0;                   0;       0];
u_max = [        150;                   1;   100e3];

% Dynamics and outputs
[dx, y] = genset_model(x, u);

%% Initial guess
% Desired engine speed
wd = rpm2rad(1500);

% Power outtake based on logistic function
P_peak = 100e3; 
T0 = 0.33; 
s = 90;
P_dem = P_peak/(1 + exp(-s*(t-T0)));

% PI Engine Speed Controller - with back calculation anti-windup
K=10; 
Ti=0.5; 
Tt=1;
I = yop.state(); % PI - integral state
e = wd - w_ice;
u_pi = K*e + I;
es = u_f - u_pi; % Non-zero when control is limited

% Simulation problem
ivp = yop.ivp( ...
    t0==0, tf==1.4, ...
    ... Controlled system
    x(t0) == x0, ... 
    der(x)== dx, ...
    u_f   == smoke_limiter(u_pi, y.u_f_max, u_min(1), u_max(1)), ...
    u_wg  == 0, ...
    P_gen == P_dem, ...
    ... Engine speed controller
    der(I) == K/Ti*e + es/Tt, ...
    I(t0)  == 0 ...
    );
sim = ivp.solve();
%%
figure(1)
subplot(411); hold on
sim.plot(t, rad2rpm(w_ice))
subplot(412); hold on
sim.plot(t, p_im)
subplot(413); hold on
sim.plot(t, p_em)
subplot(414); hold on
sim.plot(t, w_tc)

figure(2)
subplot(311); hold on
sim.plot(t, u_f)
sim.plot(t, y.u_f_max, '--', 'LineWidth', 2);
subplot(312); hold on
sim.stairs(t, u_wg)
subplot(313); hold on
sim.plot(t, P_gen)

%%
ocp = yop.ocp('Optimal Control Problem Benchmark');
ocp.min( 1e3*int(y.cylinder.fuel_massflow) );
ocp.st( ...
    ... Problem horizon
    t0==0, tf<=1.4, ...
    ... Differential constraint
    der(x) == dx, ...
    x(t0)  == x0, ...
    ... Terminal conditions
    P_gen(tf)  == 100e3, ... [W]
    int(P_gen) >= 100e3, ... [J]
    dx(tf)     == 0    , ... stationarity
    ... Box constraints
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    ... Path constraints
    y.engine.torque >= 0, ...
    hard(y.phi <= y.phi_max) ...
    );
sol = ocp.solve('intervals', 50, 'guess', sim);

%%
figure(1)
subplot(411); hold on
sol.plot(t, rad2rpm(w_ice), 'mag', 5)
subplot(412); hold on
sol.plot(t, p_im, 'mag', 5)
subplot(413); hold on
sol.plot(t, p_em, 'mag', 5)
subplot(414); hold on
sol.plot(t, w_tc, 'mag', 5)

figure(2)
subplot(311); hold on
sol.plot(t, u_f, 'mag', 5)
sol.plot(t, y.u_f_max);
subplot(312); hold on
sol.stairs(t, u_wg)
subplot(313); hold on
sol.plot(t, P_gen, 'mag', 5)

figure(3); hold on
sol.stairs(t, der(P_gen), 'mag', 5)
