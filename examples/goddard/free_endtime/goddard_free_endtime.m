%% Formulation 1
[t, t0, tf] = yop.time('t');
x = yop.state('x', 3);
u = yop.control('u');

[dx, y] = rocket_model(x, u);
rocket = y.rocket;

m_max = 215;
m_min = 68;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( rocket.height(tf) );
ocp.st(...
    ... dynamics 
    der(x) == rocket_model(x, u), ...
    ... initial condition
    t0==0, ...
    rocket.height(t0)   == 0, ...
    rocket.velocity(t0) == 0, ...
    rocket.mass(t0)     == m_max, ...
    ... Box constraints
    rocket.height   >= 0, ...
    rocket.velocity >= 0, ...
    m_min <= rocket.mass <= m_max, ...
     0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

[tt,xx,uu,pp,tx]=ocp.solve('method', 'dc', 'intervals', 70);

figure(1)
subplot(411); hold on;
plot(tx, xx(:,1))
subplot(412); hold on;
plot(tx, xx(:,2))
subplot(413); hold on;
plot(tx, xx(:,3))
subplot(414); hold on;
stairs(tt, [uu; nan])

%% Formulation 2
[t, t0, tf] = yop.time('t');
v = yop.state('v');
h = yop.state('h');
m = yop.state('m');
x = [v; h; m];

Wf = yop.control('Wf');
u = Wf;

m_max = 215;
m_min = 68;
Wf_max = 9.5;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf) );
ocp.st( ...
    t0 == 0, ...
    der(x) == rocket_model(x, u), ...    
    v(t0) == 0, ...
    h(t0) == 0, ...
    m(t0) == m_max, ...
    h >= 0, ...
    v >= 0, ...
    m_min <= m <= m_max, ...
    0 <= Wf <= Wf_max ...
    );
ocp.build().present();

dms = yop.dms(ocp, 100, 4);
sol = dms.solve();

time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(411); hold on;
plot(t_sol, x_sol(:,1))
subplot(412); hold on;
plot(t_sol, x_sol(:,2))
subplot(413); hold on;
plot(t_sol, x_sol(:,3))
subplot(414); hold on;
stairs(t_sol, [u_sol; nan])

%% Formulation 3
% Time
[t, t0, tf] = yop.time('t');

% Rocket parameters
D0 = 0.01227; beta = 0.145e-3; c = 2060;

% Constants
g0 = 9.81; r0 = 6.371e6;

% Rocket model
h  = yop.state('h');  % Rocket height
v  = yop.state('v');  % Rocket speed
m  = yop.state('m');  % Rocket mass
Wf = yop.control('Wf');  % Rocket fuel massflow

% Drag force and gravitational acceleration
F_D = D0 * exp(-beta*h) * v^2;
g   = g0*(r0/(r0+h))^2;

% Mass boundaries
m_min = 68; m_max = 215; 

% Control boundaries
Wfmin = 0; Wfmax = 9.5;

% Optimal control problem
ocp = yop.ocp();
ocp.max( h(tf) );
ocp.st( ...
    ... Rocket model
    der(h) == v              , ...
    der(v) == (Wf*c-F_D)/m-g , ...
    der(m) == -Wf            , ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_max, ...
    ... Box constraints
    h >= 0, ...
    v >=0, ...
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

ocp.build().present();

dms = yop.dms(ocp, 100, 4);
sol = dms.solve();

time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(411); hold on;
plot(t_sol, x_sol(:,1))
subplot(412); hold on;
plot(t_sol, x_sol(:,2))
subplot(413); hold on;
plot(t_sol, x_sol(:,3))
subplot(414); hold on;
stairs(t_sol, [u_sol; nan])


















