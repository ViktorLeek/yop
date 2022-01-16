# Yop
A Matlab Toolbox for Optimal Control.

#### Requirements
1. [Matlab](https://se.mathworks.com)
2. [Casadi](https://web.casadi.org/get/)

## A Few Examples
### Bryson-Denham Problem
```matlab
yops Times: t t0 tf % Parsed by position: t, t0, tf
yops States: x v    % position, speed
yops Control: a     % acceleration
yops Param: l       % maximum cart position

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) ); % int: integrate expression over problem horizon
ocp.st( t0==0, tf==1 );
ocp.st( der(v) == a );
ocp.st( der(x) == v );
ocp.st(  v(t0) == -v(tf) == 1 );
ocp.st(  x(t0) ==  x(tf) == 0 );
ocp.st( x <= l == 1/9 );
sol = ocp.solve();

figure(1);
subplot(311); hold on
sol.plot(t, x);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.stairs(t, a);
```
### Goddard's Rocket Problem
```matlab
yops Times: t t0 tf States: v h m Control: Wf

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
ocp.max( h(tf) );
ocp.st( h(t0)==0, v(t0)==0, m(t0)==m_max );
ocp.st( der(v) == (Wf*c-F_D)/m-g );
ocp.st( der(h) == v );
ocp.st( der(m) == -Wf );
ocp.st( m_min <= m  <= m_max );
ocp.st( Wfmin <= Wf <= Wfmax );

sol = ocp.solve('intervals', 50);

figure(1);
subplot(411); hold on
sol.plot(t, v);
sol.plot(t, der(h), '--'); % Test that v == der(h)
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.stairs(t, Wf);
```
### Low Thrust Orbit
```matlab
yops Times: t t0 tf
yops States: p f g h k L w scaling: [1e7,0.1,1,1,0.1,10,1]
yops Control: u size: [3,1]
yops Parameter: tau

p0 = 21837080.052835;
h0 = -0.25396764647494;
x     = [  p;  f;  g;  h;  k;     L;   w];
x_min = [2e7; -1; -1; -1; -1;    pi; 0.1];
x_max = [6e7; +1; +1; +1; +1; 18*pi;   1];
x0    = [ p0;  0;  0; h0;  0;    pi;   1];

[dx, v, Qr, r] = spacecraft(x, u, tau);

%% Initial guess
sim = yop.ivp(t0==0, tf==90e3);
sim.add( der(x) == dx );
sim.add(  x(t0) == x0 );
sim.add(  u(t)  == Qr'*v/norm(v) );
sim.add(  tau   == -25 );
res = sim.solve('solver', 'idas', 'points', 2000);

%% Optimal control problem
final = @(expr) expr(tf);

ocp = yop.ocp('Low-thrust orbit');
ocp.max( w(tf) );
ocp.st( 50e3 <= tf <= 100e3 );
ocp.st( der(x) == dx );
ocp.st(  x(t0) == x0 );
ocp.st( norm(u)^2 == 1 );
ocp.st( x_min <=  x  <= x_max );
ocp.st(   -1  <=  u  <= 1 );
ocp.st(  -50  <= tau <= 0 );
ocp.st( p(tf) == 40007346.015232 );
ocp.st( final( f^2 + g^2 ) == 0.73550320568829^2 );
ocp.st( final( h^2 + k^2 ) == 0.61761258786099^2 );
ocp.st( final( f*h + g*k ) == 0 );
ocp.st( -3 <= final( g*h - k*f ) <= 0 );

sol = ocp.solve('guess', res, 'intervals', 125, 'degree', 3);

figure(1);
subplot(321); hold on
sol.plot(t/3600, p*1e-6)
xlabel('t [h]')
ylabel('p [1e6 ft]')

subplot(322); hold on
sol.plot(t/3600, f)
xlabel('t [h]')
ylabel('f')

subplot(323); hold on
sol.plot(t/3600, g)
xlabel('t [h]')
ylabel('g')

subplot(324); hold on
sol.plot(t/3600, h)
xlabel('t [h]')
ylabel('h')

subplot(325); hold on
sol.plot(t/3600, k)
xlabel('t [h]')
ylabel('k')

subplot(326); hold on
sol.plot(t/3600, L/2/pi)
xlabel('t [h]')
ylabel('L [rad / 2\pi]')

figure(2)
subplot(311); hold on
sol.plot(t/3600, u(1))
xlabel('t [h]')
ylabel('u_r')

subplot(312); hold on
sol.plot(t/3600, u(2))
xlabel('t [h]')
ylabel('u_\theta')

subplot(313); hold on
sol.plot(t/3600, u(3))
xlabel('t [h]')
ylabel('u_h')

[ex,ey,ez] = sphere();
er = 20902000; % Earth radius [ft]

figure(3); hold on; grid on
sol.plot3(r(1), r(2), r(3), 'b', 'LineWidth', 2) 
sol.plot3(r(1).at(t0), r(2).at(t0), r(3).at(t0), 'o', 'LineWidth', 3)
sol.plot3(r(1).at(tf), r(2).at(tf), r(3).at(tf), 'o', 'LineWidth', 3)
surf(ex*er, ey*er, ez*er, 'FaceAlpha', 0.5, 'linestyle', 'none')
title('Trajectory')
xlabel('x [ft]')
ylabel('y [ft]')
zlabel('z [ft]')
```
