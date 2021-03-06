yops Times: t t0 tf
yops State: x size: [6,1] nominal: [1e3,1e3,1e3,1e2,1,1]
yops Ctrls: u size: [2,1] int: 2
yops Param: p             nominal: 0.1

x_max = [+1000; +1000; 1000; 350; +75*pi/180; +0.5*pi];
x_min = [-1000; -1000;    0;  10; -75*pi/180; -3.0*pi];

u_max = [+1.5; +75*pi/180];
u_min = [-0.5; -75*pi/180];

p_max = 0.15;
p_min = 0.005;

[dx, y] = soaring(x, u, p);

%% Optimal control problem
ocp = yop.ocp('Dynamic Soaring Problem');
ocp.min( 1e2*p );
ocp.st( ...
    t0 == 0, ...
    1 <= tf <= 30, ...
    der(x) == dx, ...
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    p_min <= p <= p_max, ...
    x(1:3).at(t0) == 0, ...
    x(1:3).at(tf) == 0, ...
    x(4:5).at(t0) == x(4:5).at(tf), ...
    x(6).at(t0) == x(6).at(tf) + 2*pi ...
    );
sol = ocp.solve('ival', 100, 'dx', 3);

%% Plot solution
figure(1); 
subplot(4,2,[1,2]); hold on
sol.plot(x(1), x(2)) 
xlabel('x [ft]')
ylabel('y [ft]')

subplot(423); hold on
sol.plot(t, x(3))
xlabel('t [s]')
ylabel('h [ft]')

subplot(424); hold on
sol.plot(t, x(4))
xlabel('t [s]')
ylabel('v [ft/s]')

subplot(425); hold on
sol.plot(t, x(5)*180/pi)
xlabel('t [s]')
ylabel('\gamma [deg]')

subplot(426); hold on
sol.plot(t, x(6)*180/pi)
xlabel('t [s]')
ylabel('\psi [deg]')

subplot(427); hold on
sol.plot(t, u(1))
xlabel('t [s]')
ylabel('C_L [-]')

subplot(428); hold on
sol.plot(t, u(2)*180/pi)
xlabel('t [s]')
ylabel('\phi [deg]')

%% Model
function [dx, y] = soaring(x, u, p)

v     = x(4);
gamma = x(5);
psi   = x(6);
c_L   = u(1);
phi   = u(2);
beta  = p;

rho  = 0.002378;
c_D0 = 0.00873;
K    = 0.045;
g    = 32.2;
m    = 5.6;
A_f  = 45.09703;

W_x     = beta * x(3);
W_x_dot = beta * v *sin(gamma); % beta + zdot
c_D     = c_D0 + K * c_L^2;
F_D     = 0.5 * rho * A_f * c_D * v^2;
F_L     = 0.5 * rho * A_f * c_L * v^2;

% Derivative
dx  = [ ...
    v * cos(gamma) * sin(psi) + W_x                                ; ...
    v * cos(gamma) * cos(psi)                                      ; ...
    v * sin(gamma)                                                 ; ...
    -F_D/m         - g*sin(gamma)  - W_x_dot*sin(psi)*cos(gamma)   ; ...
    (F_L*cos(phi)/m - g*cos(gamma) + W_x_dot*sin(psi)*sin(gamma))/v; ...
    (F_L*sin(phi)/m                - W_x_dot*cos(psi))/cos(gamma)/v ...
    ];

% Signals and constants
y.m = m;
y.g = g;
y.F_L = F_L;

end