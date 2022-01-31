yops Times: t t0 tf 
yops States: x size: [8,1] nominal: [1e2,1e5,1e5,1e4,1e2,1,10,1e5]
yops Controls: u size: [3,1] nominal: [1e2, 1, 1e5]
upshift_param;
w_ice=x(1); w_tr=x(5); u_f=u(1); u_wg=u(2); P_gen=u(3);

%% Objective
J = @(dx, y) 1e-4 * int( der( dx(5) )^2 ) + 1e-3 * int( y.cylinder.fuel_massflow );

%% First phase
[dx1, y1, c1] = coupled_gear_first(x, u);
p1 = yop.ocp('Up-shift - Phase 1');
p1.min( J(dx1, y1) );
p1.st( t0 == 0 );
p1.st( tf == 0.7 );
p1.st( der(x) == dx1 );
p1.st(  x(t0) == x0 );
p1.st( x_min <= x <= x_max );
p1.st( u_min <= u <= u_max );
p1.st( y1.engine.torque(tf) == y1.emachine.torque(tf) );
p1.st( c1{:} );

%% Second phase
[dx2, y2, c2] = decoupled(x, u);
p2 = yop.ocp('Up-shift - Phase 2');
p2.min( J(dx2, y2) );
p2.st( tf-t0 == 0.3 ); % Phase duration
p2.st( der(x) == dx2 );
p2.st(  x_min <= x <= x_max  );
p2.st(  u_min <= u <= u_max  );
p2.st( y2.engine.torque(tf) == y2.emachine.torque(tf) );
p2.st( w_ice(tf) == w_tr(tf)*i_t(2) ); % Match transmission speed
p2.st( c2{:} );

%% Third phase
[dx3, y3, c3] = coupled_gear_second(x, u);
p3 = yop.ocp('Up-shift - Phase 3');
p3.min( J(dx3, y3) );
p3.st( tf-t0 == 0.504 ); % Phase duration
p3.st( der(x) == dx3 );
p3.st( xf_min <= x(tf) <= xf_max );
p3.st(  x_min <=   x   <= x_max  );
p3.st(  u_min <=   u   <= u_max  );
p3.st( dx3(tf) == 0 );
p3.st( c3{:} );

%% Initial guess - phase 1
e1 = max(0, y1.emachine.torque - y1.engine.torque);
kp1 = 10;
ivp1 = yop.ivp();
ivp1.add( t0==0, tf==0.7 );
ivp1.add( der(x) == dx1 );
ivp1.add(  x(t0) == x0  );
ivp1.add( u_f   == kp1*e1 );
ivp1.add( u_wg  == 1 );
ivp1.add( P_gen == 0 );
sim1 = ivp1.solve();
p1.guess = sim1;

%% Initial guess - phase 2
e2  = w_ice - w_tr*i_t(2);
kp2 = 20e3;
u_pgen = min(kp2*e2, u_max(3));
u_pgen = max(u_pgen, u_min(3));
ivp2 = yop.ivp();
ivp2.add(t0 == sim1.value(tf)      );
ivp2.add(tf == sim1.value(tf) + 0.3);
ivp2.add( der(x) == dx2 );
ivp2.add(  x(t0) == sim1.value(x(tf)) );
ivp2.add( u_f   == 0 );
ivp2.add( u_wg  == 0 );
ivp2.add( P_gen == u_pgen );
sim2 = ivp2.solve();
p2.guess = sim2;

%% Initial guess - phase 3
yops State: I % PI-controller integration state
e3 = xf_max(1) - x(1);
x03 = sim2.value(x(tf));
x03(5) = x03(1)/i_t(2); % should be synchronized
ivp3 = yop.ivp();
ivp3.add(t0 == sim2.value(tf) );
ivp3.add(tf == sim2.value(tf) + 2.0);
ivp3.add( der(x) == dx3 );
ivp3.add(  x(t0) == x03 );
ivp3.add( der(I) == e3 );
ivp3.add(  I(t0) == 0 );
ivp3.add( u_f   == min(max(0, 50*e3 + 20*I), u_max(1)) );
ivp3.add( u_wg  == 0 );
ivp3.add( P_gen == -3*x(8) ); % Returned stored energy
sim3 = ivp3.solve();
p3.guess = sim3;

%% Plot initial guess
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sim1.plot(t, x(k))
    sim2.plot(t, x(k))
    sim3.plot(t, x(k))
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sim1.plot(t, u(k))
    sim2.plot(t, u(k))
    sim3.plot(t, u(k))
end

%% OCP
ocp = p1 + p2 + p3;
[sol, sol1, sol2, sol3] = ocp.solve('ival', 35, 'dx', 5);

%% Plot solution
figure(1);
for k=1:8
    subplot(4,2,k); hold on
    sol.plot(t, x(k), 'LineWidth', 2)
    sol1.plot(t, x(k))
    sol2.plot(t, x(k))
    sol3.plot(t, x(k))
end

figure(2);
for k=1:3
    subplot(3,1,k); hold on
    sol.plot(t, u(k), 'LineWidth', 2)
    sol1.plot(t, u(k))
    sol2.plot(t, u(k))
    sol3.plot(t, u(k))
end





