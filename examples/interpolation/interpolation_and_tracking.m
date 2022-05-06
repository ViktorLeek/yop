yops Times: t t0 tf 
yops States: x size: [2,1]
yops Controls: u

T = 2*pi;
tdata = 0:0.01:T;
x1_ref = interp1(tdata, sin(2.5*tdata), t);

ocp = yop.ocp('Reference tracking and interpolation');
ocp.min( 10*int( (x(1) - x1_ref)^2 ) );
ocp.st( 0 == t0, tf == T );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [0; -1] );
ocp.st( -5 <= u <= 5 );
sol = ocp.solve('ival', 100);

figure(1);
subplot(211); hold on
sol.plot(t, x(1));
sol.plot(t, x1_ref);
subplot(212); hold on
sol.stairs(t, u);