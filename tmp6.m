[t0, tf, t, x, u] = yop.ocp_variables('nx', 3, 'nu', 1);

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( 1e-5*rocket.height(tf) );
ocp.st( ...
    t0==0, tf==212, ...
    der(x) == rocket_model(x, u), ...
    rocket.height(t0)   == 0    , ...
    rocket.velocity(t0) == 0    , ...
    rocket.mass(t0)     == 215  , ...
    68 <= rocket.mass <= 215 , ...
    0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

sol = ocp.solve('intervals', 80);

figure(1);
subplot(411); hold on
sol.plot(t, x(1));
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.plot(t, u);
