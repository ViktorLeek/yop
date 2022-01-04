yopvar times: t t0 tf
yopvar state: x size: [6,1] scaling: [1e3,1e3,1e3,1e2,1,1]
yopvar  ctrl: u size: [2,1] deg: 2
yopvar param: p scaling: 0.1

x_max = [+1000; +1000; 1000; 350; +75*pi/180; +0.5*pi];
x_min = [-1000; -1000;    0;  10; -75*pi/180; -3.0*pi];

u_max = [+1.5; +75*pi/180];
u_min = [-0.5; -75*pi/180];

p_max = 0.15;
p_min = 0.005;

[dx, y] = soaring(x, u, p);

%% Initial guess
ivp = yop.simulation(t0==0, tf==30);
ivp.add( der(x) == dx );
ivp.add(  x(t0) == [0; 0; 0; 220; 0; -1] );
ivp.add( p == 0.08 );
ivp.add( u == [0.5; 0] );
sim = ivp.solve();

%% Plot initial guess
figure(1); 
subplot(4,2,[1,2]); hold on
sim.plot(x(1), x(2)) 
subplot(423); hold on
sim.plot(t, x(3))
subplot(424); hold on
sim.plot(t, x(4))
subplot(425); hold on
sim.plot(t, x(5)*180/pi)
subplot(426); hold on
sim.plot(t, x(6)*180/pi)
subplot(427); hold on
sim.plot(t, u(1))
subplot(428); hold on
sim.plot(t, u(2)*180/pi)

%% Optimal control problem
ocp = yop.ocp('Dynamic Soaring Problem');
ocp.min( 1e2*p );
ocp.st( ...
    1 <= tf <= 30, ...
    der(x) == dx, ...
    x_min <= x <= x_max, ...
    u_min <= u <= u_max, ...
    p_min <= p <= p_max, ...
    x(1:3).at(t0) == 0, ...
    x(1:3).at(tf) == 0, ...
    x(4:5).at(t0) == x(4:5).at(tf), ...
    x(6).at(t0) == x(6).at(tf) + 2*pi ...
    );
sol = ocp.solve('intervals', 100, 'degree', 3, 'guess', sim);

%% Plot solution
figure(1); 
subplot(4,2,[1,2]); hold on
sol.plot(x(1), x(2)) 
xlabel('x [ft]')
ylabel('y [ft]')

subplot(423); hold on
sol.plot(t, x(3))
xlabel('t [s]')
ylabel('h [ft]')

subplot(424); hold on
sol.plot(t, x(4))
xlabel('t [s]')
ylabel('v [ft/s]')

subplot(425); hold on
sol.plot(t, x(5)*180/pi)
xlabel('t [s]')
ylabel('\gamma [deg]')

subplot(426); hold on
sol.plot(t, x(6)*180/pi)
xlabel('t [s]')
ylabel('\psi [deg]')

subplot(427); hold on
sol.plot(t, u(1))
xlabel('t [s]')
ylabel('C_L [-]')

subplot(428); hold on
sol.plot(t, u(2)*180/pi)
xlabel('t [s]')
ylabel('\phi [deg]')