yops Times: t t0 tf
yops State: x size: [6,1] weight: [1e4,1e4,1e4,1e3,1e3,1e3]
yops Ctrl:  u size: [2,1]

Npop = 30000;
x0   = [76; 1; 36; 2; 4; 1]*Npop/120;

R = diag([25; 250]);
ocp = yop.ocp('Two-Strain Tuberculosis');
ocp.min( 1e-3*int(x(4) + x(6) + u'*R*u) );
ocp.st( t0==0, tf==5 );
ocp.st( der(x) == tuberculosis(x, u) );
ocp.st(  x(t0) == x0 );
ocp.st( 0.00 <= x <= 30e3 );
ocp.st( 0.05 <= u <= 0.95 );
ocp.st( sum(x) == Npop );

ig = yop.guess(t0, 0, tf, 5, x, x0', u, [0.95, 0.95]);

sol = ocp.solve('intervals', 100, 'guess', ig);

figure(1)
subplot(321); hold on
sol.plot(t, x(1))
subplot(322); hold on
sol.plot(t, x(2))
subplot(323); hold on
sol.plot(t, x(3))
subplot(324); hold on
sol.plot(t, x(4))
subplot(325); hold on
sol.plot(t, x(5))
subplot(326); hold on
sol.plot(t, x(6))

figure(2)
subplot(211); hold on
sol.plot(t, u(1))
subplot(212); hold on
sol.plot(t, u(2))