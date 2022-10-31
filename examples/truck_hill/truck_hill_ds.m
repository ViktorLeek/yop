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
s0 = yop.independent0();
sf = yop.independentf();
s  = yop.independent();
Ek = yop.state('scaling', 1e6);
uf = yop.control('scaling', 1e2);

M = M1;

[f, y] = truck(s, Ek, uf, M, theta);
v2ek = @(v) 0.5 * y.m * (v/3.6)^2;
v_mean = 80/3.6;

%%
tf = int(1/y.v);
t_mean = sf/v_mean; 

ocp = yop.ocp();
ocp.min(1e-6*int(uf));
ocp.st(...
    s0==0, sf==1700, ...
    der(Ek) == f, ...
    Ek(s0)==v2ek(80), ...
    Ek(sf)>=v2ek(80), ...
    v2ek(70) <= Ek <= v2ek(90), ...
    tf <= t_mean, ...
    0 <= uf <= 300 ...
    );

sol = ocp.solve('ival', 400, 'dx', 2);

figure(1);
subplot(211); hold on
sol.plot(s, y.v*3.6)
subplot(212); hold on
sol.stairs(s, uf)


%%
function [dEk, y] = truck(s, Ek, uf, M, theta)
g   = 9.81;  % Gravitational acc [m/s^2]
m   = 40e3;  % Mass of yruck [kg]
Af  = 10;    % Frontal area [m^2]
rho = 1.292; % Density of air [kg/m^3]
cd  = 0.5;   % Drag coefficient [-]
cr  = 0.006; % Rolling res. param [-]
ig  = 3;     % Total gear ratio [-]
rw  = 0.5;   % Wheel diameter [m]
ncyl = 6;    % Number of cylinders
qlhv = 42.9e6;

v2 = 2*Ek/m;
v = sqrt(v2);

we = v / rw * ig; % Engine speed [rad/s]
Ne = we * 30/pi;  % Engine speed [rpm]
Wf = uf * Ne * ncyl * (1e-6/120); % Fuel flow [kg/s]
Me = M(uf);       % Engine torque [Nm]
Pe = pi/30*Ne*Me; % Engine power [W]
eta= Pe/Wf/qlhv;  % Engine efficiency

Fp  = Me * ig / rw; % Engine torque [Nm] -> Propulsive force [N]
Fr  = m * g * cr * cos(theta(s)); % Rolling resistance [N]
Fa  = 0.5 * rho * Af * cd * v2; % Aerodynamic resistance [N]
Fg  = m * g * sin(theta(s)); % Force due to gravity [N] 
dEk = Fp - Fr - Fa - Fg; % acceleration [N]

y.v = v;
y.m = m;
y.Ne = Ne;
y.Me = Me;
y.Pe = Pe/v;
y.Wf = Wf/v;
y.eta = eta;
y.m = m;
end