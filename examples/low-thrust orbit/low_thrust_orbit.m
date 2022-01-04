yopvar times: t t0 tf
yopvar states: p f g h k L w scaling: [1e7,0.1,1,1,0.1,10,1]
yopvar control: u size: [3,1]
yopvar param: tau

p0 = 21837080.052835;
h0 = -0.25396764647494;
x     = [  p;  f;  g;  h;  k;     L;   w];
x_min = [2e7; -1; -1; -1; -1;    pi; 0.1];
x_max = [6e7; +1; +1; +1; +1; 18*pi;   1];
x0    = [ p0;  0;  0; h0;  0;    pi;   1];

[dx, v, Qr, r] = spacecraft(x, u, tau);

%% Initial guess
ivp = yop.simulation( ...
    t0==0, tf==90e3, ...
    der(x) == dx, ...
    x(t0) == x0, ...
    u == Qr'*v/norm(v), ...
    tau == -25 ...
    );
sim = ivp.solve('solver', 'idas', 'points', 2000);

%% Plot guess
figure(1);
subplot(321); hold on
sim.plot(t, p)
subplot(322); hold on
sim.plot(t, f)
subplot(323); hold on
sim.plot(t, g)
subplot(324); hold on
sim.plot(t, h)
subplot(325); hold on
sim.plot(t, k)
subplot(326); hold on
sim.plot(t, L)

figure(2)
subplot(311); hold on
sim.plot(t, u(1))
subplot(312); hold on
sim.plot(t, u(2))
subplot(313); hold on
sim.plot(t, u(3))

%% Optimal control problem
final = @(expr) expr(tf);

ocp = yop.ocp('Low-thrust orbit');
ocp.max( w(tf) );
ocp.st( 50e3 <= tf <= 100e3 );
ocp.st( der(x) == spacecraft(x, u, tau) );
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

sol = ocp.solve('guess', sim, 'intervals', 125, 'degree', 3);

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
sol.plot3(rv(1,1)  , rv(2,1)  , rv(3,1)  , 'o', 'LineWidth', 3)
sol.plot3(rv(1,end), rv(2,end), rv(3,end), 'o', 'LineWidth', 3)
surf(ex*er, ey*er, ez*er, 'FaceAlpha', 0.5, 'linestyle', 'none')
title('Trajectory')
xlabel('x [ft]')
ylabel('y [ft]')
zlabel('z [ft]')
