%% Optimal Control Problem Benchmark
t0    = yop.time0('t0');
tf    = yop.timef('tf');
t     = yop.time('t');
w_ice = yop.state('w_ice'); % Engine angular velocity
p_im  = yop.state('p_im');  % Intake manifold pressure
p_em  = yop.state('p_em');  % Exhause manifold pressure
w_tc  = yop.state('w_tc');  % Turbocharger angular velocity
x = [w_ice; p_im; p_em; w_tc]; % State vector

u_f   = yop.control('u_f');   % Fuel injection per cycle per cylinder
u_wg  = yop.control('u_wg');  % Wastegate control 0-close, 1-fully open
P_gen = yop.control('P_gen'); % Generator power
u = [u_f; u_wg; P_gen]; % Control vector

[dx, y] = genset_model(x, u);
x0 = [rpm2rad(800); 1.0143e5; 1.0975e5; 2.0502e3];
x_min = [rpm2rad(800);  8.1e4; 9.1e4;  500];
x_max = [rpm2rad(2500); 3.5e5; 4.0e5; 15e3];
u_min = [0; 0; 0];
u_max = [150; 1; 100e3];

%% Initial guess
I    = yop.state(); % PI - integral state

K=10; Ti=1; Tt=1;
wd = rpm2rad(1500);
e = wd - w_ice;
u_pi = K*e + I;
es = u_f - u_pi;

P_peak = 120e3; T0 = 0.75; s = 6.25;
P_dem = P_peak/(1 + exp(-s*(t-T0)));

ivp = yop.ivp( ...
    t0==0, tf==3, ...
    ... Controlled system
    x(t0) == x0, ... 
    der(x) == dx, ...
    u_wg == 0, ...
    ... Smoke limiter
    u_f == smoke_limiter(u_pi, y.cylinder.fuel_max, 0, 150), ...
    ... Engine speed controller
    der(I) == K/Ti*e + es/Tt, ...
     I(t0) == 0, ...
    ... Power outtake
    P_gen == P_dem ...
    );
sol = ivp.solve();

figure(1)
subplot(411); hold on
sol.plot(t, rad2rpm(w_ice))
subplot(412); hold on
sol.plot(t, p_im)
sol.plot(t, p_em)
subplot(413); hold on
sol.plot(t, w_tc)
subplot(414); hold on
sol.plot(t, I)

figure(2)
subplot(211); hold on
sol.plot(t, u_f)
sol.plot(t, y.cylinder.fuel_max, '--', 'LineWidth', 2);
subplot(212); hold on
sol.plot(t, P_gen)

%%

ocp = yop.ocp('Optimal Control Problem Benchmark');
ocp.min( int(y.cylinder.fuel_massflow) ); % Minimize total fuel mass
ocp.st( ...
    ... Problem horizon
    t0 == 0, ...
    1.1 <= tf <= 1.4, ...
    ... Differential constraint
    der(x) == dx, ...
     x(t0) == x0, ...
    ... Terminal conditions
     P_gen(tf) == 100e3, ... [W]
    int(P_gen) == 100e3, ... [J]
    dx(1:4).at(tf) == 0    , ... stationarity
    ... Box constraints
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    ... Path constraints
    y.turbine.BSR_min <= y.turbine.BSR <= y.turbine.BSR_max, ...
    y.engine.power <= y.engine.power_limit, ...
    y.compressor.pressure_ratio <= y.compressor.surge_line, ...
    y.cylinder.fuel_to_air_ratio <= 1/y.cylinder.lambda_min ...
    );
sol = ocp.solve('intervals', 100);


