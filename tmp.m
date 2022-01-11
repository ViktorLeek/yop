yops times: t t0 tf
yops states: x size: [3, 1] weight: [1e3,1e5,1e2]
yops controls: u weight: 10

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( y.rocket.height(tf)*1e-5 );
ocp.st( t0==0 );
ocp.st( der(x) == dx );
ocp.st( y.rocket.height(t0)   == 0 );
ocp.st( y.rocket.velocity(t0) == 0 );
ocp.st( y.rocket.mass(t0)     == 215 );
ocp.st( 68 <= y.rocket.mass <= 215 );
ocp.st( 0 <= y.rocket.fuel_mass_flow <= 9.5 );
%% Initial guess: maximum thrust until fuel is out
% 
f = @(t,x) rocket_model(x, yop.IF(x(3)>=68, 9.5, 0));
[t_sim, x_sim] = ode15s(f, [0, 200], [0;0;215]);

% Reconstruct control input
u_sim = x_sim(:,3);
u_sim(u_sim < 68)  = 0;
u_sim(u_sim >= 68) = 9.5;

% Initial guess, pairwise (variable, value)
guess = yop.guess(t, t_sim, x, x_sim, u, u_sim);

figure(1);
subplot(411); hold on
plot(t_sim, x_sim(:,1))
subplot(412); hold on
plot(t_sim, x_sim(:,2))
subplot(413); hold on
plot(t_sim, x_sim(:,3))
subplot(414); hold on
plot(t_sim, u_sim);

%% Solve OCP

opts = struct;
opts.ipopt.print_level = 0;
sol = ocp.solve('intervals', 50, 'guess', guess, 'opts', opts);

figure(1);
subplot(411); hold on
sol.plot(t, x(1), 'mag', 5);
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.stairs(t, u);
