%% Bang-bang Control
%   http://www.ee.ic.ac.uk/ICLOCS/ExampleBangBang2.html
yops Times: t t0 tf
yops States: x size: [2,1]
yops Control: u

% Dynamics: dx = Ax + Bu
A = [0,1;0,0];
B = [0; 1];

% Bounds phase 1
u_min = -2; u_max =  1;
x_max = [300; +200];
x_min = [-10; -200];
xf    = [300;    0];

% There
p1 = yop.ocp('Bang-bang - Phase 1');
p1.min( 0 );
p1.st( 0 == t0 <= tf );
p1.st( der(x) == A*x + B*u );
p1.st(  x(t0) ==  0 );
p1.st(  x(tf) == xf );
p1.st(   x_min <= x  <=  x_max );
p1.st(   u_min <= u  <=  u_max );

% Bounds phase 2
x_max = [310; +200];
x_min = [  0; -200];
u_max =  2;
u_min = -1; 

% and back again
p2 = yop.ocp('Bang-bang - Phase 2');
p2.min( tf );
p2.st( t0 <= tf );
p2.st( der(x) == A*x + B*u );
p2.st(  x(tf) ==  0 );
p2.st(   x_min <= x  <=  x_max );
p2.st(   u_min <= u  <=  u_max );

% A hobbit's tale
ocp = p1 + p2;
[sol, sol1, sol2] = ocp.solve('intervals', 15, 'degree', 2);

% by Bilbo Baggins
figure(1); 
subplot(311); hold on;
sol.plot(t, x(1))
sol2.plot(t, x(1))
subplot(312); hold on;
sol.plot(t, x(2))
sol2.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)
sol2.stairs(t, u)
