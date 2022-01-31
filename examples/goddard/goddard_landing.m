%% Goddard's Rocket Problem - Including landing phase
yops Times: t t0 tf
yops States: v h m nominal: [1e3,1e5,1e2]
yops Ctrls: Wf int: 2 nominal: 10

% Parameters
D0 = 0.01227; beta = 0.145e-3; c = 2060;

% Constants
g0 = 9.81; r0 = 6.371e6;

% Drag force and gravitational acceleration
F_D = sign(v) * D0 * exp(-beta*h) * v^2;
g   = g0*(r0/(r0+h))^2;

% Mass boundaries
m_min = 68; m_max = 215; 

% Control boundaries
Wfmin = 0; Wfmax = 9.5;

% Optimal control problem
p1 = yop.ocp();
p1.max( h(tf)*1e-4 );
p1.st( t0==0 );
p1.st( h(t0)==0, v(t0)==0, m(t0)==m_max );
p1.st( der(v) == (Wf*c-F_D)/m-g );
p1.st( der(h) == v );
p1.st( der(m) == -Wf );
p1.st( m_min <= m  <= m_max );
p1.st( Wfmin <= Wf <= Wfmax );
sol1 = p1.solve(); % Using solution to first phase as initial guess for the phase
p1.guess = sol1; 

p2 = yop.ocp();
p2.min( 0 );
p2.st( tf >= t0 );
p2.st( der(v) == (Wf*c-F_D)/m-g );
p2.st( der(h) == v );
p2.st( der(m) == -Wf );
p2.st(  v(tf) == 0 );
p2.st(  h(tf) == 0 );
p2.st( h >= 0 );
p2.st( m_min <= m  <= m_max );
p2.st( Wfmin <= Wf <= Wfmax );

ocp = p1 + p2;

sol = ocp.solve('ival', 100, 'dx', 5);

figure(1);
subplot(411); hold on
sol.plot(t, v, 'mag', 4);
subplot(412); hold on
sol.plot(t, h, 'mag', 4);
subplot(413); hold on
sol.plot(t, m, 'mag', 4);
subplot(414); hold on
sol.plot(t, Wf, 'mag', 4);
