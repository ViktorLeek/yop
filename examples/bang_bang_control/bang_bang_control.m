%% Bang-bang Control
%   From: http://www.ee.ic.ac.uk/ICLOCS/ExampleBangBang.html
yops Times: t t0 tf
yops States: x size: [2,1]
yops Control: u

A = [0,1;0,0];
B = [0; 1];
tf_min = 0; tf_max = 35;
u_min = -2; u_max =  1;
x_max = [300; +200];
x_min = [-10; -200];
xf    = [300;    0];

ocp = yop.ocp('Bang-bang');
ocp.min( 0.1*tf );
ocp.st( t0 == 0 );
ocp.st( der(x) == A*x + B*u );
ocp.st(  x(t0) ==  0 );
ocp.st(  x(tf) == xf );
ocp.st(  tf_min <= tf <= tf_max );
ocp.st(   x_min <= x  <=  x_max );
ocp.st(   u_min <= u  <=  u_max );
sol = ocp.solve('intervals', 30, 'degree', 2);

figure(1); 
subplot(311); hold on;
sol.plot(t, x(1))
subplot(312); hold on;
sol.plot(t, x(2))
subplot(313); hold on;
sol.stairs(t, u)
