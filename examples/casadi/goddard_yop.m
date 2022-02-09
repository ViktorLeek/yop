yops Times: t t0 tf
yops States: x size: [3,1] 
yops Controls: u

[dx, y] = rocket_model(x, u);

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( y.rocket.height(tf) );
ocp.st( t0==0 );
ocp.st( der(x) == dx );
ocp.st( y.rocket.height(t0)   == 0 );
ocp.st( y.rocket.velocity(t0) == 0 );
ocp.st( y.rocket.mass(t0)     == 215 );
ocp.st( 68 <= y.rocket.mass <= 215 );
ocp.st( 0 <= y.rocket.fuel_mass_flow <= 9.5 );
sol = ocp.solve();

figure(1)
subplot(311); hold on
sol.plot(t, x(1));
subplot(312); hold on
sol.plot(t, x(2));
subplot(313); hold on
sol.plot(t, x(3));

figure(2); hold on
sol.stairs(t, u);



















