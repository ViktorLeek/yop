%% Optimal Control Problem Benchmark

[t, t0, tf] = yop.time();

w_ice = yop.state('w_ice'); % Engine angular velocity
p_im  = yop.state('p_im');  % Intake manifold pressure
p_em  = yop.state('p_em');  % Exhause manifold pressure
w_tc  = yop.state('w_tc');  % Turbocharger angular velocity
E_gen = yop.state('E_gen'); % Generator energy
x = [w_ice; p_im; p_em; w_tc; E_gen]; % State vector

u_f   = yop.control('u_f');   % Fuel injection per cycle per cylinder
u_wg  = yop.control('u_wg');  % Wastegate control 0-close, 1-fully open
P_gen = yop.control('P_gen'); % Generator power
u = [u_f; u_wg; P_gen]; % Control vector

[dx, y] = genset_model(x, u);

dx14 = dx(1:4);

ocp = yop.ocp('Optimal Control Problem Benchmark');
ocp.min( int(y.cylinder.fuel_massflow) ); % Minimize total fuel mass
ocp.st( ...
    ... Problem horizon
    t0 == 0, ...
    1.1 <= tf <= 1.4, ...
    ... Differential constraint
    der(x) == dx, ...
    ... Initial conditions
    w_ice(t0) == rpm2rad(800),... [rad/s]
     p_im(t0) == 1.0143e5    ,... [Pa]
     p_em(t0) == 1.0975e5    ,... [Pa]
     w_tc(t0) == 2.0502e3    ,... [rad/s]
    E_gen(t0) == 0           ,... [J]
    ... Terminal conditions
    E_gen(tf) == 100e3, ... [J]
    P_gen(tf) == 100e3, ... [W]
     dx14(tf) == 0    , ... stationarity
    ... Box constraints
    rpm2rad(800) <= w_ice <= rpm2rad(2500), ...
    8.1e4 <= p_im  <= 3.5e5 , ...
    9.1e4 <= p_em  <= 4.0e5 , ...
      500 <= w_tc  <= 15e3  , ...
        0 <= E_gen <= 3.0e6 , ...
        0 <= u_f   <= 150   , ...
        0 <= u_wg  <= 1     , ...
        0 <= P_gen <= 100e3 , ...
    ... Path constraints
    y.turbine.BSR_min <= y.turbine.BSR <= y.turbine.BSR_max, ...
    y.engine.power <= y.engine.power_limit, ...
    y.compressor.pressure_ratio <= y.compressor.surge_line, ...
    y.cylinder.fuel_to_air_ratio <= 1/y.cylinder.lambda_min ...
    );
ocp.build().present();
dms = yop.dms(ocp, 25, 4);
sol = dms.solve();

%% Optimal Control Problem Benchmark
w_ice = yop.state('w_ice'); % Engine angular velocity
p_im  = yop.state('p_im');  % Intake manifold pressure
p_em  = yop.state('p_em');  % Exhause manifold pressure
w_tc  = yop.state('w_tc');  % Turbocharger angular velocity
E_gen = yop.state('E_gen'); % Generator energy
x = [w_ice; p_im; p_em; w_tc; E_gen]; % State vector

u_f   = yop.control('u_f');   % Fuel injection per cycle per cylinder
u_wg  = yop.control('u_wg');  % Wastegate control 0-close, 1-fully open
P_gen = yop.control('P_gen'); % Generator power
u = [u_f; u_wg; P_gen]; % Control vector

w_ice.m_value = casadi.MX.sym('v');
p_im.m_value = casadi.MX.sym('v');
p_em.m_value = casadi.MX.sym('v');
w_tc.m_value = casadi.MX.sym('v');
E_gen.m_value = casadi.MX.sym('v');

u_f.m_value = casadi.MX.sym('v');
u_wg.m_value = casadi.MX.sym('v');
P_gen.m_value = casadi.MX.sym('v');

dx = genset_model(x,u);

xc = x.evaluate;
uc = u.evaluate;
dx = dx.evaluate;

w = [xc(:); uc(:)];
j = casadi.Function('j', {w}, {jacobian(dx, w)});
j([rpm2rad(800); 1.0143e5; 1.0975e5; 2.0502e3; 0; 50; 0; 0])

%%
x = casadi.MX.sym('x', 5, 1);
u = casadi.MX.sym('u', 3, 1);
dx = genset_model(x,u);
w = [x(:); u(:)];
j = casadi.Function('j', {w}, {jacobian(dx, w)});
j([rpm2rad(800); 1.0143e5; 1.0975e5; 2.0502e3; 0; 50; 0; 0])


