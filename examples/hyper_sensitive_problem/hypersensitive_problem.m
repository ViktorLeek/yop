% Hypersensitive problem from GPOPS-II
yops Times: t t0 tf State: x Control: u

ocp = yop.ocp('Hyper-senstive problem');
ocp.min( 0.5*int(x^2 + u^2) );
ocp.st( t0==0, tf==1000, der(x) == -x^3 + u, x(t0)==1.5, x(tf)==1 );
sol = ocp.solve('ival', 500);

figure(1)
subplot(211); hold on
sol.plot(t, x, 'mag', 5)
subplot(212); hold on
sol.stairs(t, u)