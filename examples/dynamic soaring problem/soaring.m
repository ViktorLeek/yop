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