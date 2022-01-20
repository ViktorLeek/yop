%% Multp-phase Bryson-Denham Problem
yops Times: t t0 tf % Parsed by position: t, t0, tf
yops States: x v    % position, speed
yops Controls: a       % acceleration

p1 = yop.ocp('Bryson - Phase 1');
p1.min( 1/2 * int(a^2) );
p1.st( t0 == 0 );
p1.st( tf == 0.5 );
p1.st( der(x) == v );
p1.st( der(v) == a );
p1.st( x(t0)  == 0 );
p1.st( v(t0)  == 1 );
p1.st( x <= 1/9 );

p2 = yop.ocp('Bryson - Phase 2');
p2.min( 1/2 * int(a^2) );
p2.st( tf == 1 );
p2.st( der(x) == v );
p2.st( der(v) == a );
p2.st( x(tf)  == 0 );
p2.st( v(tf)  == -1 );
p2.st( x <= 1/9 );

ocp = p1 + p2;
[sol, s1, s2] = ocp.solve('intervals', [5, 5], 'degree', [4, 4]);

figure(1);
subplot(311); hold on
sol.plot(t, x, 'LineWidth', 2);
s1.plot(t, x);
s2.plot(t, x);
subplot(312); hold on
sol.plot(t, v, 'LineWidth', 2);
s1.plot(t, v);
s2.plot(t, v);
subplot(313); hold on
% sol.stairs(t, a, 'LineWidth', 2);
% s1.stairs(t, a);
% s2.stairs(t, a);
% sol.plot(t, a, 'LineWidth', 2);
s1.plot(t, a, 'x-');
s2.plot(t, a, 'x-');