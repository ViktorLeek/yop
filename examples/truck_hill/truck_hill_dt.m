%% Data
uf =    [254 232 205 193 166 140 114 101 88 75 62 49 37 25]';
Me = 10*[239 219 194 182 156 130 104  91 78 65 52 39 26 13]';
A1 = [0*uf+1 uf ];
A2 = [0*uf+1 uf  uf.^2];
B  = Me;
c1 = A1\B;
c2 = A2\B;

%% Torque model
M1 = @(uf)               c1(2)*uf + c1(1);
M2 = @(uf) c2(3)*uf.^2 + c2(2)*uf + c2(1);

%% Road slope
theta0 = atan(15/300);
L = @(s, s0) theta0/(1 + exp(-0.1*(s-s0)));
theta  = @(s) L(s, 400) - L(s, 700) - L(s, 1100) + L(s, 1400);

%% Variables
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
s  = yop.state();
v  = yop.state();
u  = yop.control();
x  = [s; v];

%% Guess

% Mean velocity and constraints
v_des = 80/3.6; % [m/s]
v_max = 90/3.6; % [m/s]
v_min = 70/3.6; % [m/s]

% Contorl input limits
u_max = 300; % [mg/cyc/cyl]
u_min = 0;   % [mg/cyc/cyl]

% Total distance
d = 1700;        % [m]
t_max = d/v_des; % [s]

% Dynamics and signals
[f, y] = truck(x, u, M1, theta);

%% Simulation
% Controller parameters
e  = v_des - v; % [m/s]
K  = 1000;      % Proportional constant
Ti = 1;         % Integral time constant
Ts = 0.1;       % Tracking time constant

PI = PI_controller(K, Ti, Ts, u_min, u_max);

% Simulation
ivp = yop.ivp(t0==0, tf==t_max); % Simulation time
ivp.add( der(x) == f );          % 
ivp.add(  x(t0) == [0; v_des] );
ivp.add( PI.u == u  );
ivp.add( PI.e == v_des - v );
ivp.add( PI.dynamics );

% Simulate
sim = ivp.solve('points', 150, 'solver', 'idas');

%%
figure(1);
subplot(211); hold on
sim.plot(s, v*3.6)
subplot(212); hold on
sim.plot(s, u)

%% OCP - M1
ocp = yop.ocp();
ocp.min( int(y.Wf) );
ocp.st(0 == t0 <= tf <= t_max);
ocp.st( der(x) == f );
ocp.st( s(t0) == 0 );
ocp.st( s(tf) == d );
ocp.st( v(t0) == v_des );
ocp.st( v(tf) >= v_des );
ocp.st( v_min <= v <= v_max);
ocp.st( u_min <= u <= u_max);
sol1 = ocp.solve('ival', 400, 'guess', sim);

%% OCP - M2

[f, y] = truck(x, u, M2, theta);
ocp = yop.ocp();
ocp.min( int(y.Wf) );
ocp.st(0 == t0 <= tf <= t_max);
ocp.st( der(x) == f );
ocp.st( s(t0) == 0 );
ocp.st( s(tf) == d );
ocp.st( v(t0) == v_des );
ocp.st( v(tf) >= v_des );
ocp.st( v_min <= v <= v_max);
ocp.st( u_min <= u <= u_max);
sol2 = ocp.solve('ival', 400, 'guess', sim);


%%
figure(1);
subplot(211); hold on
sol1.plot(s, v*3.6)
sol2.plot(s, v*3.6)

subplot(212); hold on
sol1.stairs(s, u)
sol2.stairs(s, u)

%%
function [dx, y] = truck(x, uf, M, theta)
g    = 9.81;  % Gravitational acc [m/s^2]
m    = 40e3;  % Mass of yruck [kg]
Af   = 10;    % Frontal area [m^2]
rho  = 1.292; % Density of air [kg/m^3]
cd   = 0.5;   % Drag coefficient [-]
cr   = 0.006; % Rolling res. param [-]
ig   = 3;     % Total gear ratio [-]
rw   = 0.5;   % Wheel diameter [m]
nr   = 2;     % Revolutions per cycle [-]
ncyl = 6;     % Number of cylinders [-]

s = x(1); % Position [m]
v = x(2); % Velocity [m/s]

Ne  = v * ig / rw / 2 / pi; % Engine speed [rps]
Wf  = uf * Ne * ncyl / nr * 1e-6; % Fuel flow [kg/s]
Me  = M(uf); % [Nm]

Fp = Me * ig / rw; % Engine torque [Nm] -> Propulsive force [N]
Fr = m * g * cr * cos(theta(s)); % Rolling resistance [N]
Fa = 0.5 * rho * Af * cd * v^2; % Aerodynamic resistance [N]
Fg = m * g * sin(theta(s));
dv = (Fp - Fr - Fa - Fg)/m;
dx = [v; dv];
y.Wf = Wf; % y - struct for signals
end

function pi = PI_controller(K, Ti, Tt, u_min, u_max)
e  = yop.algebraic(); % PI input
u  = yop.algebraic(); % PI output
I  = yop.state();     % integration state
uu = K*e + I;         % Unsaturated control
us = min(max(uu, u_min), u_max); % Upper bound
es = us-uu;           % Saturation error
pi.u = u; 
pi.e = e;
pi.dynamics = [u; der(I); I(yop.time0)] == [us; K/Ti*e + es/Tt; 0];
end