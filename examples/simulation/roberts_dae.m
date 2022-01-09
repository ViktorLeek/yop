%% DAE Simulation - Roberts problem
%   The Robertson problem found in hb1ode.m
%   https://se.mathworks.com/help/matlab/ref/ode15s.html
yops times: t t0 tf states: x1 x2 algebraic: z

sim = yop.ivp(t0==4e-6, tf==4e6);
sim.add( der(x1) == -0.04 * x1 + 1e4 * x2 * z );
sim.add( der(x2) == +0.04 * x1 - 1e4 * x2 * z - 3e7 * x2^2 );
sim.add(    0    == x1 + x2 + z - 1 );
sim.add( x1(t0)  == 1 );
sim.add( x2(t0)  == 0 );
sim.add(  z(t0)  == 0 );
sol_15s = sim.solve('solver', 'ode15s', 'reltol', 1e-4, 'abstol', [1e-6, 1e-10, 1e-6]);
sol_idas = sim.solve('solver', 'idas', 'opts', struct('grid', 4*logspace(-6,6)));

figure(1);
sol_15s.semilogx(t, x1); hold on
sol_15s.semilogx(t, 1e4*x2)
sol_15s.semilogx(t, z)

set(gca,'ColorOrderIndex',1)
sol_idas.semilogx(t, x1, 'x-'); hold on
sol_idas.semilogx(t, 1e4*x2, 'x-');
sol_idas.semilogx(t, z, 'x-');