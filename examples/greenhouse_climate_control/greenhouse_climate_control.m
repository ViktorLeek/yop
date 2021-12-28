[t0, tf, t, x, u] = yop.vars('nx', 2, 'nu', 1);

p = [7.5e-8; 1; 0.1; 4.55e-4; 136.4];
I = max(0, 800*sin(4*pi*t/tf - 0.65*pi));
T0 = 15 + 10*sin(4*pi*t/tf - 0.65*pi);

ocp = yop.ocp('Greenhouse Climate Control');
ocp.min( -p(5)*x(1).at(tf) + int(p(4)*u) );
ocp.st( t0==0, tf==48 );
ocp.st( der(x(1)) == p(1)*I*x(2) );
ocp.st( der(x(2)) == p(2)*(T0-x(2)) + p(3)*u );
ocp.st( x(t0) == [0; 10] );
ocp.st( 0 <= u <= 10 );
sol = ocp.solve('intervals', 20);

figure(1)
subplot(311); hold on;
sol.plot(t, 1200*x(1))
subplot(312); hold on;
sol.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)