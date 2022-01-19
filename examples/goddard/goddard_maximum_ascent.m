yops Times: t t0 tf
yops States: v h m
yops Ctrls: T

% Parameters
c=0.5; g0=1; h0=1; D0=310; b=500;

% Drag
g = g0*(h0/h)^2;
D = D0*exp(-b*h);
F_D = D*v^2;

% Optimal control problem
ocp = yop.ocp();
ocp.max( h(tf)^2 );
ocp.st( t0==0, h(t0)==1, v(t0)==0, m(t0)==1 );
ocp.st( der(v) == (T - F_D)/m - g );
ocp.st( der(h) == v );
ocp.st( der(m) == -T/c );
ocp.st( h >= 1, m >= 0.6 );
ocp.st( 0 <= T <= 3.5 );

sol = ocp.solve('intervals', 50);

figure(1);
subplot(411); hold on
sol.plot(t, v);
subplot(412); hold on
sol.plot(t, h);
subplot(413); hold on
sol.plot(t, m);
subplot(414); hold on
sol.plot(t, T);