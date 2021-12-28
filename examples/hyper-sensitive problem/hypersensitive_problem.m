% Hypersensitive problem from GPOPS-II
yopvar t0 tf t x u

ocp = yop.ocp('Hyper-senstive problem');
ocp.min( 0.5*int(x^2 + u^2) );
ocp.st( der(x) == -x^3 + u, x(t0)==1.5, x(tf)==1, tf==1000 );
sol = ocp.solve('intervals', 500);

figure(1)
subplot(211); hold on
sol.plot(t, x, 'mag', 5)
subplot(212); hold on
sol.stairs(t, u)