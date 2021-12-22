%% ODE with parameter
%   https://se.mathworks.com/help/matlab/ref/ode15s.html
yop_time t t0 tf
yop_parameter A B
y = yop.state(2);

sim = yop.simulation();
sim.add(t0==0, tf==5);
sim.add(der(y) == [y(2); A/B*t*y(1)]);
sim.add( y(t0) == [0; 0.01]);
sim.add(A==1, B==2);
sol = sim.solve('solver', 'cvodes');
sol.plot(t, y(1), '-o', t, y(2), '-.')