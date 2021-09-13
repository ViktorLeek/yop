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

%%
load('data_MVEM2.mat')
load('benchmark_data.mat')
[t, t0, tf] = yop.time();

w_ice = yop.state('w_ice');
p_im  = yop.state('p_im');
p_em  = yop.state('p_em');
w_tc  = yop.state('w_tc');
x     = [w_ice; p_im; p_em; w_tc];

u_f   = yop.control('u_f');
u_wg  = yop.control('u_wg');
P_gen = yop.control('P_gen');
u     = [u_f; u_wg; P_gen];

[dx, c] = MVEM2(w_ice, p_im, p_em, w_tc, u_f, u_wg, P_gen, param);

ocp = yop.ocp('Optimal Control Benchmark');
ocp.min( (int(c.dot_m_f))^2 );
ocp.st( ...
    t0 == 0, ...
    1.1 <= tf <= 1.4, ...
    der(x) == dx, ...
    x(t0) == x0', ...
    param.x_min <= x <= param.x_max, ...
    param.u_min <= u <= param.u_max, ...
    c.P_GEN <= 100000/param.control_norm(3), ...
    c.P_GEN(tf) >= 100000/param.control_norm(3), ...
    int(c.P_GEN) >= 100000/param.control_norm(3), ...
    c.P_ice <= param.P_ice_max, ...
    c.P_ice <= param.cPice(1)*power(c.W_ICE,2)+param.cPice(2)*c.W_ICE+param.cPice(3), ...
    c.P_ice <= param.cPice(4)*power(c.W_ICE,2)+param.cPice(5)*c.W_ICE+param.cPice(6), ...
    0 <= c.phi <= 1/param.lambda_min, ...
    param.BSR_min <= c.BSR <= param.BSR_max, ...
    c.Pi_c <= param.c_mc_surge(1)*c.dot_m_c_corr+param.c_mc_surge(2) ...
    );
ocp.build();

%%

N = 50;
d = 4;
dc = yop.dc(ocp, N, d, 'legendre');
nlp = dc.build();


tt = [];
for n=1:N
    for r=1:d+1
        tt = [tt(:); dc.t{n,r}];
    end
end
tt = [tt(:); dc.t{N+1,1}];
time = casadi.Function('x', {dc.t0, dc.tf}, {tt});
t00 = 0;
tf0 = 1;
x_grid = full(time(t00, tf0));
x_guess = interp1(int_ex.states(:,1), int_ex.states(:,2:end), x_grid');
x_guess = x_guess';
x00 = x_guess(:);
u_grid = linspace(t00, tf0, dc.N);
u_guess = interp1(int_ex.control(:,1), int_ex.control(:,2:end), u_grid');
u_guess = u_guess';
u00 = u_guess(:);
w0 = [t00; tf0; x00; u00];

prob = struct; prob.f = nlp.f; prob.x = nlp.x; prob.g = nlp.g;
solver = casadi.nlpsol('solver', 'ipopt', prob);
sol = solver( ...
    'x0', w0, ...
    'lbx', nlp.x_lb, ...
    'ubx', nlp.x_ub, ...
    'ubg', nlp.g_ub, ...
    'lbg', nlp.g_lb ...
    );

%%
tt = [];
for n=1:N
    for r=1:d+1
        tt = [tt(:); dc.t{n,r}];
    end
end
tt = [tt(:); dc.t{N+1,1}];
time = casadi.Function('x', {dc.w}, {tt});
tx_sol = full(time(sol.x));

tt = [];
for n=1:N+1
    tt = [tt(:); dc.t{n,1}];
end
time = casadi.Function('x', {dc.w}, {tt});
tu_sol = full(time(sol.x));

xx = [];
for xk=dc.x
    xx = [xx; xk.y'];
end
state = casadi.Function('x', {dc.w}, {xx});
x_sol = full(state(sol.x));

control = casadi.Function('x', {dc.w}, {horzcat(dc.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(411); hold on;
plot(tx_sol, x_sol(:,1).*param.state_norm(1))
subplot(412); hold on;
plot(tx_sol, x_sol(:,2).*param.state_norm(2))
subplot(413); hold on;
plot(tx_sol, x_sol(:,3).*param.state_norm(3))
subplot(414); hold on;
plot(tx_sol, x_sol(:,4).*param.state_norm(4))

figure(2)
subplot(311);
stairs(tu_sol, [u_sol(:,1); nan].*param.control_norm(1))
subplot(312);
stairs(tu_sol, [u_sol(:,2); nan].*param.control_norm(2))
subplot(313);
stairs(tu_sol, [u_sol(:,3); nan].*param.control_norm(3))


