%% Constraint rocket acceleration below a certain height
hh = 0:1e1:2e5;
h0 = 0.2e5;
k = 1e-2;
c = @(h) 1-1./(1+exp(-k*(h-h0)));

plot(hh,c(hh))
%%
yops Times: t t0 tf nominal: [1e2,1e0,1e2]
yops States: v h m  nominal: [1e3,1e5,1e2]
yops Controls: u    nominal: 10 int: 1

x     = [  v;   h;   m];
x0    = [  0;   0; 215];
x_max = [inf; inf; 215];
x_min = [  0;   0;  68];

u_max = 9.5;
u_min = 0.0;

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf)*1e-5 );
ocp.st( der(x) == dx, ...
         x(t0) == x0, t0==0, ...
        x_min <= x <= x_max, ...
        u_min <= u <= u_max );
ocp.hard( der(v)*c(h) <= 20 );

sol = ocp.solve();
%%
figure(1);
subplot(411); hold on
sol.plot(t, v);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.plot(t, u, 'mag', 5);

figure(2); hold on;
sol.plot(t, der(v)*c(h))
sol.plot(t, der(v), '--')

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
y.rocket.acceleration   = dv;
y.rocket.velocity       = v;
y.rocket.height         = h;
y.rocket.mass           = m;
y.rocket.fuel_mass_flow = F;
y.drag.coefficient      = D;
y.drag.force            = F_D;
y.gravity               = g;
end