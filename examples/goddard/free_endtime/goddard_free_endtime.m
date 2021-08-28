%% Formulation 1
[t, t0, tf] = yop.independent('t');
x = yop.state('x', 3);
u = yop.control('u');

[dx, y] = rocket_model(x, u);
rocket = y.rocket;

m0 = 215;
mf = 68;

ocp = yop.ocp();
ocp.max(y.rocket.height(tf));
ocp.st(...
    ... dynamics 
    der(x) == rocket_model(x, u), ...
    ... initial condition
    t0  == 0, ...
    rocket.height(t0)   == 0, ...
    rocket.velocity(t0) == 0, ...
    rocket.mass(t0)     == m0, ...
    ... Box constraints
    rocket.height   >= 0, ...
    rocket.velocity >= 0, ...
    mf <= rocket.mass <= m0, ...
     0 <= rocket.fuel_mass_flow <= 9.5 ...
    );
ocp.build();
ocp.present();

%% Formulation 2
[t, t0, tf] = yop.time('t');
v = yop.state('v');
h = yop.state('h');
m = yop.state('m');
F = yop.control('F');

m_max = 215;
m_min = 68;
F_max = 9.5;

[dx, y] = rocket_model([v; h; m], F);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf) );
ocp.st( ...
    t0 == 0, ...
    der([v; h; m]) == dx, ...
    v(t0) == 0, ...
    h(t0) == 0, ...
    m(t0) == m_max, ...
    h >= 0, ...
    v >= 0, ...
    m_min <= m <= m_max, ...
    0 <= F <= F_max ...
    );
ocp.build();
ocp.present();

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
grp = yop.ocp();
grp.max( h(tf) );
grp.st( ...
    ... Rocket model
    der(h) == v              , ...
    der(v) == (Wf*c-F_D)/m-g , ...
    der(m) == -Wf            , ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_min, ...
    ... Box constraints
    h >= 0, ...
    v >=0, ...
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

grp.build();
grp.present();






















