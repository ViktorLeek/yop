%% Formulation 1
[t0, tf, t, x, u] = yop.ocp_variables('nx', 3, 'nu', 1);

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
% ocp.max( rocket.height(tf) );
ocp.max( int(rocket.velocity) );
ocp.st( ...
    der(x) == rocket_model(x, u), ...
    rocket.height(t0)   == 0    , ...
    rocket.velocity(t0) == 0    , ...
    rocket.mass(t0)     == 215  , ...
    68 <= rocket.mass <= 215, ...
    0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

sol = ocp.solve('intervals', 40);

figure(1);
subplot(411); hold on
sol.plot(t, x(1));
subplot(412); hold on
sol.plot(t, x(2));
sol.plot(t(t==42), x(2).at(t==42), 'x')
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.stairs(t, u);

% sol.save('Goddard.mat');
% sol = yop.load('Goddard', t, t0, tf, x, u, p);

%% Formulation 2
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');
v  = yop.state('v');
h  = yop.state('h');
m  = yop.state('m');
Wf = yop.control('Wf');

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
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');

% Rocket model
h  = yop.state('h');  % Rocket height
v  = yop.state('v');  % Rocket speed
m  = yop.state('m');  % Rocket mass
Wf = yop.control('Wf');  % Rocket fuel massflow

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

% Optimal control problem
ocp = yop.ocp();
ocp.max( h(tf) );
ocp.st( ...
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


















