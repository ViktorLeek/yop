%% ODE with parameter
%   https://se.mathworks.com/help/matlab/ref/ode15s.html
yopvar times: t t0 tf state: y size: [2, 1] parameters: A B

sim = yop.ivp();
sim.add(t0==0, tf==5);
sim.add(der(y) == [y(2); A/B*t*y(1)]);
sim.add( y(t0) == [0; 0.01]);
sim.add(A==1, B==2);
sol = sim.solve('solver', 'cvodes');
sol.plot(t, y(1), '-o', t, y(2), '-.')