yops Times: t t0 tf
yops States: r u v
yops Control: theta

% Parameters
T    = 0.1405;
m0   = 1;
mdot = 0.07489;

% Expressions
a = T/(m0 - abs(mdot)*t);
G = [u; v-1/sqrt(r)];
d2r = pi/180;

% Optimal control problem
ocp = yop.ocp('Max Radius Orbit Transfer');
ocp.max( r(tf) );
ocp.st( t0==0, tf==3.3155 );
ocp.st( der(r) == u );
ocp.st( der(u) == v^2/r - 1/r^2 + a*sin(theta) );
ocp.st( der(v) == -u*v/r + a*cos(theta) );
ocp.st(  r(t0) == 1 );
ocp.st(  u(t0) == 0 );
ocp.st(  v(t0) == 1 );
ocp.st(  G(tf) == 0  );
sol = ocp.solve();

% Project control on the ival [0, 2*pi]
theta_norm = mod(theta+1000*2*pi, 2*pi);

% Present results
disp(['Maximum final radius: ', num2str(sol.value(r(tf)))])

figure(1);
subplot(511); hold on
sol.plot(t, r);
subplot(512); hold on
sol.plot(t, u);
subplot(513); hold on
sol.plot(t, v);
subplot(514); hold on
sol.plot(t, theta_norm);
subplot(515); hold on
sol.plot(t, G)
