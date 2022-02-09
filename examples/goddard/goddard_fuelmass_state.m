yops Times: t t0 tf
yops States: x size: [3,1] 
yops Control: u

[dx, rocket] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( rocket.height(tf) );
ocp.st( t0==0 );
ocp.st( der(x) == dx );
ocp.st( rocket.height(t0) == 0 );
ocp.st( rocket.speed(t0) == 0 );
ocp.st( rocket.fuel_mass(t0) == 150 );
ocp.st( rocket.fuel_mass >= 0 );
ocp.st( 0 <= rocket.fuel_massflow <= 9.5 );
sol = ocp.solve();

figure(1)
subplot(411); hold on
sol.plot(t, rocket.height); title('height')
subplot(412); hold on
sol.plot(t, rocket.speed); title('speed')
subplot(413); hold on
sol.plot(t, rocket.fuel_mass); title('fuel mass')
subplot(414); hold on
sol.plot(t, u); title('control'); xlabel('time')


function [dx, y] = rocket_model(x, u)
% States and control
h = x(1); v = x(2); mf = x(3); F  = u;

% Parameters
m0   = 68;
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
m = m0 + mf;
dh = v;
dv = (F*c-F_D)/m-g;
dmf = -F;
dx = [dh; dv; dmf];

% Signals y
y = struct('height',h,'speed',v,'fuel_mass',mf,'fuel_massflow',F);
end 
