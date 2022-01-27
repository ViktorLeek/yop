%% Implementation 1
yops Times: t t0 tf
yops States: x size: [3,1] nominal: [1e3,1e5,1e2]
yops Controls: u nominal: 10

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( y.rocket.height(tf)*1e-5 );
ocp.st( t0==0 );
ocp.st( der(x) == dx );
ocp.st( y.rocket.height(t0)   == 0 );
ocp.st( y.rocket.velocity(t0) == 0 );
ocp.st( y.rocket.mass(t0)     == 215 );
ocp.st( 68 <= y.rocket.mass <= 215 );
ocp.st( 0 <= y.rocket.fuel_mass_flow <= 9.5 );

sol = ocp.solve();

figure(1);
subplot(411); hold on
sol.plot(t, x(1));
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.plot(t, u);

%% Implementation 2 
% PW quadratic control (by integration), integration of velocity for height
yops Times: t t0 
yops State: x1 x2 x3 weight: [1e3,1e5,1e2] 
yops Ctrls: u int: 2
x = [x1; x2; x3];

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( int(rocket.velocity)*1e-5 );
ocp.st( t0==0 );
ocp.st( der(x) == rocket_model(x, u) );
ocp.st( rocket.height(t0) == 0 );
ocp.st( rocket.velocity(t0) == 0 );
ocp.st( rocket.mass(t0) == 215 );
ocp.st( 68 <= rocket.mass <= 215 );
ocp.st(  0 <= rocket.fuel_mass_flow <= 9.5 );

sol = ocp.solve('ival', 50);

figure(1);
subplot(411); hold on
sol.plot(t, x(1));
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.plot(t, u, 'mag', 5); % Exact upsampling 5 times

%% Implementation 2 - piecewise polynomials control input
yops Times: t t0 tf weight: [1e2,1e0,1e2]
yops States: v h m  weight: [1e3,1e5,1e2]
yops Controls: u    weight: 10

x     = [  v;   h;   m];
x0    = [  0;   0; 215];
x_max = [inf; inf; 215];
x_min = [  0;   0;  68];

u_max = 9.5;
u_min = 0.0;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf)*1e-5 );
ocp.st( t0==0 );
ocp.st( der(x) == rocket_model(x, u) );
ocp.st(  x(t0) == x0 );
ocp.st( x_min <= x <= x_max );
ocp.st( u_min <= u <= u_max );

% Alternative way of entering options.
% 100 control intervals, state polynomial degree 5, control polynomial
% degree 2, control trajectory is continuous (not its derivative),
% maximum number of iterations 1000, acceptable tolerance 1e-7
yoptions Ival: 100 State: 5 Control: 2 -control_continuity
yoptions Maxiter: 1000 Acctol: 1e-7
sol = ocp.solve(yoptions);

figure(1);
subplot(411); hold on
sol.plot(t, v);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.plot(t, u, 'mag', 5);

figure(2); % derivatives
subplot(311); hold on
sol.plot(t, u, 'mag', 5);
subplot(312); hold on
sol.plot(t, der(u), 'mag', 5);
subplot(313); hold on
sol.plot(t, der(der(u)), 'mag', 5);

%% Implementation 3 - Formulate dynamics into the problem directly
yops Times: t t0 tf
yops States: v h m scaling: [1e3,1e5,1e2]
yops Ctrls: Wf int: 2 scaling: 10

% Parameters
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
ocp.max( h(tf)*1e-5 );
ocp.st( t0==0 );
ocp.st( h(t0)==0, v(t0)==0, m(t0)==m_max );
ocp.st( der(v) == (Wf*c-F_D)/m-g );
ocp.st( der(h) == v );
ocp.st( der(m) == -Wf );
ocp.st( m_min <= m  <= m_max );
ocp.st( Wfmin <= Wf <= Wfmax );

sol = ocp.solve('ival', 50);

figure(1);
subplot(411); hold on
sol.plot(t, v);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.plot(t, Wf);

%% Model

function [dx, y] = rocket_model(x, u)
% States and control
v = x(1);
h = x(2);
m = x(3);
F = u;

% Parameters
D0   = 0.01227;
beta = 0.145e-3;
c    = 2060;
g0   = 9.81;
r0   = 6.371e6;

% Drag and gravity
D   = D0*exp(-beta*h);
F_D = D*v^2;
g   = g0*(r0/(r0+h))^2;

% Dynamics
dv = (F*c-F_D)/m-g;
dh = v;
dm = -F;
dx = [dv;dh;dm];

% Signals y
y.rocket.velocity       = v;
y.rocket.height         = h;
y.rocket.mass           = m;
y.rocket.fuel_mass_flow = F;
y.drag.coefficient      = D;
y.drag.force            = F_D;
y.gravity               = g;
end 

