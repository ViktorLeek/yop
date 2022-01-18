%% ODE Simulation - Van der Pol equations for mu=1 
%   The van der pol model from vdp1.m
%   https://se.mathworks.com/help/matlab/ref/ode45.html
yops Times: t t0 tf State: x size: [2,1]
sim = yop.ivp(t0==0, tf==20);
sim.add(der(x) == [x(2); (1-x(1)^2)*x(2)-x(1)]);
sim.add(x(t0)  == [2; 0]);
% res = sim.solve('solver', 'idas', 'points', 200);
sol = sim.solve('solver', 'ode15s');

figure(1); hold on;
sol.plot(t, x(1));
sol.plot(t, x(2));