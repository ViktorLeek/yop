%% Third Order State Constraints
%   From: https://www.bocop.org/robbins/
yops Times: t t0 tf 
yops States: y1 y2 y3
yops Control: u
y = [y1; y2; y3];

a = 3;
ocp = yop.ocp('Robbins');
ocp.min( 1/2 * int( a*y1 + 1/2*u^2 ) );
ocp.st( t0==0, tf==10 );
ocp.st( der(y) == [y2;y3;u] );
ocp.st(  y(t0) == [ 1;-2;0] );
ocp.st(  y1 >= 0 );
sol = ocp.solve('intervals', 100, 'degree', 3);

figure(1); 
subplot(411); hold on;
sol.plot(t, y1)
subplot(412); hold on;
sol.plot(t, y2)
subplot(413); hold on;
sol.plot(t, y3)
subplot(414); hold on;
sol.stairs(t, u)