yops Times: t t0 tf State: x Control: u

ocp = yop.ocp('Isopermetric Constraint');
ocp.min( int(x) );
ocp.st( t0==0, tf==1 );
ocp.st( x(t0)==1, x(tf)==0 );
ocp.st( der(x) == -sin(x) + u );
ocp.st( -10 <= x <= 10 );
ocp.st(  -4 <= u <=  4 );
ocp.st( int(u^2) == 10 );
sol = ocp.solve('ival', 25);

figure(1)
subplot(211); hold on;
sol.plot(t, x)
subplot(212); hold on;
sol.stairs(t, u)