yopvar times: t t0 tf
yopvar state: x        size: [6,1] scaling: [1e3,1e3,1e3,1e2,1,1]
yopvar  ctrl: u        size: [2,1]
yopvar param: p                    scaling: 0.1

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
ivp.add( u == [0.3 + 0.5*sin(0.1*t); (-40 + 10*sin(0.55*t))*pi/180] );
sim = ivp.solve('solver', 'ode15s');

%% Plot path
sim.plot3(x(1), x(2), x(3))

%%
ocp = yop.ocp('Dynamic Soaring Problem');
ocp.min( p );
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
sol = ocp.solve('intervals', 50, 'degree', 5, 'guess', sim);

%%
% sol.plot3(x(1), x(2), x(3))
% sol.plot(t, u(2))
sol.plot(t, x(6))
