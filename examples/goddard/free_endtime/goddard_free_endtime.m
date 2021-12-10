%% Formulation 1
[t0, tf, t, x, u] = yop.ocp_variables('nx', 3, 'nu', 1);

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( y.rocket.height(tf) );
ocp.st( ...
    t0==0, tf==212, ...
    der(x) == dx, ...
    y.rocket.height(t0)   == 0    , ...
    y.rocket.velocity(t0) == 0    , ...
    y.rocket.mass(t0)     == 215  , ...
    y.rocket.velocity(t<50) <= 500, ...
    68 <= y.rocket.mass <= 215 , ...
    0 <= y.rocket.fuel_mass_flow <= 9.5 ...
    );

sol = ocp.solve('intervals', 50);

figure(1);
subplot(411); hold on
sol.plot(t, x(1), 'mag', 5);
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.plot(t, u);

% sol.save('Goddard.mat');
% sol = yop.load('Goddard', t, t0, tf, x, u, p);

%% Formulation 1 - variation 1: PWX control, integration of velocity
[t0, tf, t, x] = yop.ocp_variables('nx', 3);
u = yop.control('pw', 'quadratic'); % <-- Option for control parametrization

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( int(rocket.velocity) );
ocp.st( ...
    der(x) == rocket_model(x, u), ...
    rocket.height(t0)   == 0    , ...
    rocket.velocity(t0) == 0    , ...
    rocket.mass(t0)     == 215  , ...
    68 <= rocket.mass <= 215 , ...
    0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

sol = ocp.solve('intervals', 80);

figure(1);
subplot(411); hold on
sol.plot(t, x(1), 'mag', 5);
subplot(412); hold on
sol.plot(t, x(2), 'mag', 5);
subplot(413); hold on
sol.plot(t, x(3), 'mag', 5);
subplot(414); hold on
sol.plot(t, u, 'mag', 5);

%% Formulation 2
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
v  = yop.state();
h  = yop.state();
m  = yop.state();
Wf = yop.control();

x = [v; h; m];
u = Wf;

m_max = 215;
m_min = 68;
Wf_max = 9.5;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf) );
ocp.st( ...
    der(x) == rocket_model(x, u), ...    
    v(t0) == 0, ...
    h(t0) == 0, ...
    m(t0) == m_max, ...
    m_min <= m <= m_max, ...
    0 <= Wf <= Wf_max ...
    );

sol = ocp.solve('intervals', 50);

figure(1);
subplot(411); hold on
sol.plot(t, v, 'mag', 2);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.stairs(t, Wf);

%% Formulation 3
% Time
t0 = yop.t0();
tf = yop.tf();
t  = yop.t();

% States
h  = yop.state();  % Rocket height
v  = yop.state();  % Rocket speed
m  = yop.state();  % Rocket mass
Wf = yop.control();  % Rocket fuel massflow

% Rocket parameters
D0 = 0.01227; beta = 0.145e-3; c = 2060;

% Constants
g0 = 9.81; r0 = 6.371e6;

% Drag force and gravitational acceleration
F_D = D0 * exp(-beta*h) * v^2;
g   = g0*(r0/(r0+h))^2;

% Mass boundaries
m_min = 68; m_max = 215; 

% Control boundaries
Wfmin = 0; Wfmax = 9.5;

% dh = der(h);
% M(2:4) = f(x)

% Optimal control problem
ocp = yop.ocp();
ocp.max( h );
ocp.st( ...
    ... tf == 212, ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_max, ...
    ... Rocket model
    der(v) == (Wf*c-F_D)/m-g , ...
    der(h) == v              , ...
    der(m) == -Wf            , ...
    ... Box constraints
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

sol = ocp.solve('intervals', 50);

figure(1);
subplot(411); hold on
sol.plot(t, v);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.stairs(t, Wf);


















