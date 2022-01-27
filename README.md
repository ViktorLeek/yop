# Yop
A Matlab Toolbox for Optimal Control.

#### Requirements
1. [Matlab](https://se.mathworks.com)
2. [Casadi](https://web.casadi.org/get/)

## A Few Examples
### Bryson-Denham Problem
```matlab
yops Times: t t0 tf
yops States: x size: [2,1]
yops Controls: u

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(u^2) );
ocp.st( 0 == t0 < tf == 1 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve();

figure(1);
subplot(211); hold on
sol.plot(t, x);
subplot(212); hold on
sol.stairs(t, u);
```
### Goddard's Rocket Problem
```matlab
yops Times: t t0 tf
yops States: x size: [3,1]
yops Controls: u

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( y.rocket.height(tf) );
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
sol.stairs(t, u);
```
### Low Thrust Orbit
```matlab
yops Times: t t0 tf
yops States: p f g h k L w nominal: [1e7,0.1,1,1,0.1,10,1]
yops Control: u size: [3,1]
yops Param: tau

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

%% Plot guess
figure(1);
subplot(321); hold on
res.plot(t, p)
subplot(322); hold on
res.plot(t, f)
subplot(323); hold on
res.plot(t, g)
subplot(324); hold on
res.plot(t, h)
subplot(325); hold on
res.plot(t, k)
subplot(326); hold on
res.plot(t, L)

figure(2)
subplot(311); hold on
res.plot(t, u(1))
subplot(312); hold on
res.plot(t, u(2))
subplot(313); hold on
res.plot(t, u(3))

%% Optimal control problem
final = @(expr) expr(tf);

ocp = yop.ocp('Low-thrust orbit');
ocp.max( w(tf) );
ocp.st( t0 == 0 );
ocp.st( 50e3 <= tf <= 100e3 );
ocp.st( der(x) == dx );
ocp.st(  x(t0) == x0 );
ocp.st(   u'*u == 1 );
ocp.st( x_min <=  x  <= x_max );
ocp.st(   -1  <=  u  <= 1 );
ocp.st(  -50  <= tau <= 0 );
ocp.st( p(tf) == 40007346.015232 );
ocp.st( final( f^2 + g^2 ) == 0.73550320568829^2 );
ocp.st( final( h^2 + k^2 ) == 0.61761258786099^2 );
ocp.st( final( f*h + g*k ) == 0 );
ocp.st( -3 <= final( g*h - k*f ) <= 0 );

sol = ocp.solve('guess', res, 'ival', 125, 'dx', 3);

%% Plot solution
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

rv = sol.value(r);
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

function [dx, v_vec, Qr, r_vec] = spacecraft(x, u, tau)

p=x(1); f=x(2); g=x(3); h=x(4); k=x(5); L=x(6); w=x(7);
ur=u(1); ut=u(2); uh=u(3);

T   = 4.446618e-3;    % [lb]
Isp = 450;            % [s]
mu  = 1.407645794e16; % [ft^3/s^2]
gs  = 32.174;         % [ft/s^2]
Re  = 20925662.73;    % [ft]
J2  = 1082.639e-6;
J3  = -2.565e-6;
J4  = -1.608e-6;

% Gravitational disturbing acceleration
q   = 1 + f*cos(L) + g*sin(L);
r   = p/q;
a2  = h^2 - k^2;
chi = sqrt(h^2 + k^2);
s2  = 1 + chi^2;

rx = (r/s2) * (cos(L) + a2*cos(L) + 2*h*k*sin(L));
ry = (r/s2) * (sin(L) - a2*sin(L) + 2*h*k*cos(L));
rz = (2*r/s2) * (h*sin(L) - k*cos(L));
r_vec = [rx; ry; rz];
r_mag = norm(r_vec);

vx = -(1/s2) * sqrt(mu/p) * (  sin(L) + a2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + a2*g);
vy = -(1/s2) * sqrt(mu/p) * ( -cos(L) + a2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + a2*f);
vz =  (2/s2) * sqrt(mu/p) * (h*cos(L) +  k*sin(L) + f*h + g*k);
v_vec = [vx; vy; vz];

r_Xv = cross(r_vec, v_vec);
r_Xv_mag = norm(r_Xv);
r_Xv_Xr  = cross(r_Xv, r_vec);

ir1 = rx/r_mag;
ir2 = ry/r_mag;
ir3 = rz/r_mag;
ir = [ir1; ir2; ir3];
it1 = r_Xv_Xr(1) / (r_Xv_mag * r_mag);
it2 = r_Xv_Xr(2) / (r_Xv_mag * r_mag);
it3 = r_Xv_Xr(3) / (r_Xv_mag * r_mag);
it = [it1; it2; it3];
ih1 = r_Xv(1) / r_Xv_mag;
ih2 = r_Xv(2) / r_Xv_mag;
ih3 = r_Xv(3) / r_Xv_mag;
ih = [ih1; ih2; ih3];
enir = ir3;

enirir1 = enir * ir1;
enirir2 = enir * ir2;
enirir3 = enir * ir3;
enenirir1 = 0 - enirir1;
enenirir2 = 0 - enirir2;
enenirir3 = 1 - enirir3;
enenirir_mag = norm([enenirir1, enenirir2, enenirir3]);
in1 = enenirir1 / enenirir_mag;
in2 = enenirir2 / enenirir_mag;
in3 = enenirir3 / enenirir_mag;

% Geocentric latitude
sinphi = rz/norm([rx, rz]);
cosphi = sqrt(1 - sinphi^2);

% Legendre polynomials
P2  = 0.5000 * ( 3 * sinphi^2 - 2);
P3  = 0.5000 * ( 5 * sinphi^3 - 3*sinphi);
P4  = 0.1250 * (35 * sinphi^4 - 30*sinphi^2 + 3);
dP2 = 3.0000 * sinphi;
dP3 = 0.5000 * ( 15 * sinphi - 3);
dP4 = 0.1250 * (140 * sinphi^3 - 60*sinphi);


% Oblate earth perturbations
sumn =   (Re/r)^2 * dP2*J2 +   (Re/r)^3 * dP3*J3 +   (Re/r)^4 * dP4*J4;
sumr = 3*(Re/r)^2 *  P2*J2 + 4*(Re/r)^3 *  P3*J3 + 5*(Re/r)^4 *  P4*J4;
Dgn = -mu*cosphi/(r^2) * sumn;
Dgr = -mu/(r^2)*sumr;
Dgnin1 = Dgn * in1;
Dgnin2 = Dgn * in2;
Dgnin3 = Dgn * in3;
Dgrir1 = Dgr * ir1;
Dgrir2 = Dgr * ir2;
Dgrir3 = Dgr * ir3;
Dg1 = Dgnin1 - Dgrir1;
Dg2 = Dgnin2 - Dgrir2;
Dg3 = Dgnin3 - Dgrir3;
Deltag1 = ir(1)*Dg1 + ir(2)*Dg2 + ir(3)*Dg3;
Deltag2 = it(1)*Dg1 + it(2)*Dg2 + it(3)*Dg3;
Deltag3 = ih(1)*Dg1 + ih(2)*Dg2 + ih(3)*Dg3;

% Thrust acceleration
DeltaT1 = gs*T*(1 + 0.01*tau)/w * ur;
DeltaT2 = gs*T*(1 + 0.01*tau)/w * ut;
DeltaT3 = gs*T*(1 + 0.01*tau)/w * uh;

% Total acceleration
Delta1 = Deltag1 + DeltaT1;
Delta2 = Deltag2 + DeltaT2;
Delta3 = Deltag3 + DeltaT3;

% Differential equations of motion
dp = 2*p/q * sqrt(p/mu) * Delta2;
df =  sqrt(p/mu) * sin(L)*Delta1 + sqrt(p/mu)*(1/q)*((q+1)*cos(L) + f)*Delta2 - sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L))*Delta3;
dg = -sqrt(p/mu) * cos(L)*Delta1 + sqrt(p/mu)*(1/q)*((q+1)*sin(L) + g)*Delta2 + sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L))*Delta3;
dh =  sqrt(p/mu) * (s2*cos(L)/(2*q))*Delta3;
dk =  sqrt(p/mu) * (s2*sin(L)/(2*q))*Delta3;
dL =  sqrt(p/mu) * (1/q)*(h*sin(L) - k*cos(L))*Delta3 + sqrt(mu*p)*((q/p)^2);
dw = -(T*(1+0.01*tau)/Isp);
dx = [dp; df; dg; dh; dk; dL; dw];
Qr = [ir, it, ih];
```

### A multi-phase problem with simulated initial guess - Min jerk in driveline
```matlab
yops Times: t t0 tf
yops States: x size: [8,1] nominal: [1e2,1e5,1e5,1e4,1e2,1,10,1e5]
yops Controls: u size: [3,1] nominal: [1e2, 1, 1e5]
upshift_param;
w_ice=x(1); w_tr=x(5); u_f=u(1); u_wg=u(2); P_gen=u(3);

%% Objective
J = @(dx, y) 1e-4 * int( der( dx(5) )^2 ) + 1e-3 * int( y.cylinder.fuel_massflow );

%% First phase
[dx1, y1, c1] = coupled_gear_first(x, u);
p1 = yop.ocp('Up-shift - Phase 1');
p1.min( J(dx1, y1) );
p1.st( t0 == 0 );
p1.st( tf == 0.7 );
p1.st( der(x) == dx1 );
p1.st(  x(t0) == x0 );
p1.st( x_min <= x <= x_max );
p1.st( u_min <= u <= u_max );
p1.st( y1.engine.torque(tf) - y1.emachine.torque(tf) == 0 );
p1.st( c1{:} );

%% Second phase
[dx2, y2, c2] = decoupled(x, u);
p2 = yop.ocp('Up-shift - Phase 2');
p2.min( J(dx2, y2) );
p2.st( tf-t0 == 0.3 ); % Phase duration
p2.st( der(x) == dx2 );
p2.st(  x_min <= x <= x_max  );
p2.st(  u_min <= u <= u_max  );
p2.st( y2.engine.torque(tf) - y2.emachine.torque(tf) == 0 );
p2.st( w_ice(tf) == w_tr(tf)*i_t(2) ); % Match transmission speed
p2.st( c2{:} );

%% Third phase
[dx3, y3, c3] = coupled_gear_second(x, u);
p3 = yop.ocp('Up-shift - Phase 3');
p3.min( J(dx3, y3) );
p3.st( tf-t0 == 0.504 ); % Phase duration
p3.st( der(x) == dx3 );
p3.st( xf_min <= x(tf) <= xf_max );
p3.st(  x_min <=   x   <= x_max  );
p3.st(  u_min <=   u   <= u_max  );
p3.st( dx3(tf) == 0 );
p3.st( c3{:} );

%% Initial guess - phase 1
e1 = max(0, y1.emachine.torque - y1.engine.torque);
kp1 = 10;
ivp1 = yop.ivp();
ivp1.add( t0==0, tf==0.7 );
ivp1.add( der(x) == dx1 );
ivp1.add(  x(t0) == x0  );
ivp1.add( u_f   == kp1*e1 );
ivp1.add( u_wg  == 1 );
ivp1.add( P_gen == 0 );
sim1 = ivp1.solve();
p1.guess = sim1;

%% Initial guess - phase 2
e2  = w_ice - w_tr*i_t(2);
kp2 = 20e3;
u_pgen = min(kp2*e2, u_max(3));
u_pgen = max(u_pgen, u_min(3));
ivp2 = yop.ivp();
ivp2.add(t0 == sim1.value(tf)      );
ivp2.add(tf == sim1.value(tf) + 0.3);
ivp2.add( der(x) == dx2 );
ivp2.add(  x(t0) == sim1.value(x(tf)) );
ivp2.add( u_f   == 0 );
ivp2.add( u_wg  == 0 );
ivp2.add( P_gen == u_pgen );
sim2 = ivp2.solve();
p2.guess = sim2;

%% Initial guess - phase 3
yops State: I % PI-controller integration state
e3 = xf_max(1) - x(1);
x03 = sim2.value(x(tf));
x03(5) = x03(1)/i_t(2); % should be synchronized
ivp3 = yop.ivp();
ivp3.add(t0 == sim2.value(tf) );
ivp3.add(tf == sim2.value(tf) + 2.0);
ivp3.add( der(x) == dx3 );
ivp3.add(  x(t0) == x03 );
ivp3.add( der(I) == e3 );
ivp3.add(  I(t0) == 0 );
ivp3.add( u_f   == min(max(0, 50*e3 + 20*I), u_max(1)) );
ivp3.add( u_wg  == 0 );
ivp3.add( P_gen == -3*x(8) ); % Returned stored energy
sim3 = ivp3.solve();
p3.guess = sim3;

%% Plot initial guess
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sim1.plot(t, x(k))
    sim2.plot(t, x(k))
    sim3.plot(t, x(k))
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sim1.plot(t, u(k))
    sim2.plot(t, u(k))
    sim3.plot(t, u(k))
end

%% OCP
ocp = p1 + p2 + p3;
[sol, sol1, sol2, sol3] = ocp.solve('ival', 15, 'dx', 5);

%% Plot solution
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sol.plot(t, x(k), 'mag', 5, '--', 'LineWidth', 2)
    sol1.plot(t, x(k), 'mag', 5)
    sol2.plot(t, x(k), 'mag', 5)
    sol3.plot(t, x(k), 'mag', 5)
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sol.plot(t, u(k), '--', 'LineWidth', 2)
    sol1.plot(t, u(k))
    sol2.plot(t, u(k))
    sol3.plot(t, u(k))
end

```
