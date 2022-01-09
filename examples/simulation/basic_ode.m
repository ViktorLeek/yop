%% Basic ode
yops times: t t0 tf state: x
sim = yop.ivp(t0==0, tf==2);
sim.add( der(x) == -10*t );
sim.add(  x(t0) == 1 ) ;
cv_sol = sim.solve('solver', 'idas');
ode_sol = sim.solve('solver', 'ode45');
figure(1); hold on;
cv_sol.plot(t, x);
ode_sol.plot(t, x);

%% Short representation
%   Using default solver IDAS
yop.ivp(t0==0, tf==2, der(x)==-10*t, x(t0)==1).solve().plot(t, x);