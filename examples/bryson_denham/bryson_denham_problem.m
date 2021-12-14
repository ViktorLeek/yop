%% Original formulation
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state();   % position
v  = yop.state();   % speed
a  = yop.control(); % acceleration
l  = yop.parameter(); % maximum cart position

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( t0==0, tf==1 );
ocp.st( der(v) == a );
ocp.st( der(x) == v );
ocp.st( v(t0) == -v(tf) == 1 );
ocp.st( x(t0) ==  x(tf) == 0 );
ocp.st( x <= l == 1/9 );
sol = ocp.solve();

figure(1);
subplot(311); hold on
sol.plot(t, x);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.plot(t, a);

%% State vector, minimum value and traveled distance
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state(2); % [position, speed]
u  = yop.control('pw', 'quadratic'); % acceleration

J = 1/2 * int(u^2);
ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( J );
ocp.st(   t0==0, tf==1    );
ocp.st( der(x) == [x(2); u]  );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve('intervals', 50);

J_min = sol.value(J);
dist = sol.value(int(abs(x(2))));

disp(['Minimum cost is ', num2str(J_min)]);
disp(['Traveled distance is ', num2str(dist)]);

figure(1);
subplot(211); hold on
sol.plot(t, x);
subplot(212); hold on
sol.plot(t, u);

%% Guaranteed box constraints for boundary conditions
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state();     % position
v  = yop.state();     % speed
a  = yop.control();   % acceleration

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == v, ...
    der(v) == a, ...
    x(t0) == 0, ...
    v(t0) == 1, ...
    x(tf) == 0, ...
    v(tf) == -1, ...
    x <= 1/9 ...
    );

sol = ocp.solve('intervals', 15, 'degree', 2);
figure(1);
subplot(311); hold on
sol.plot(t, x, 'mag', 5);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.stairs(t, a);

%% Compact representation
[t0, tf, t, x, u] = yop.ocp_variables('nx', 2, 'nu', 1);
ocp = yop.ocp().min( 1/2 * int(u^2) );
ocp.st(   t0==0, tf==1    );
ocp.st( der(x) == [x(2); u]  );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve('intervals', 20, 'degree', 2);
figure(1); hold on
sol.plot(t, x);
sol.stairs(t, u);

%% Minreal
[t0, tf, t, x, u] = yop.ocp_variables('nx', 2, 'nu', 1);
yop.ocp().min(1/2*int(u^2)).st(tf==1, der(x)==[x(2);u], x(t0)==[0; 1], ...
    x(tf)==[0;-1], x(1)<=1/9).solve().plot(t, [x;u]);

%% Trade-off between control effort and traveled distance
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state(); % position
v  = yop.state(); % speed
a  = yop.control(); % acceleration
l  = yop.parameter(); % maximum position of the cart

beta = 30; % trade-off factor

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) + beta*l );
ocp.st( t0==0, tf==1 );
ocp.st( der(v) == a );
ocp.st( der(x) == v );
ocp.st( v(t0) == -v(tf) == 1 );
ocp.st( x(t0) == x(tf) == 0 );
ocp.st( x <= l );

sol = ocp.solve('intervals', 20, 'degree', 2);

figure(1);
subplot(311); hold on
sol.plot(t, x, 'mag', 5);
sol.plot(t, 0*t + l);
subplot(312); hold on
sol.plot(t, v, 'mag', 5);
subplot(313); hold on
sol.stairs(t, a, 'mag', 5);

%% Simulation
t0 = yop.time0();
tf = yop.timef();
t  = yop.time();
x  = yop.state();     % position
v  = yop.state();     % speed
a  = yop.control();   % acceleration
l  = yop.parameter(); % maximum position of the cart

ivp = yop.ivp( ...
    t0==0, tf==1, ...
    x(t0) == 0  , ...
    v(t0) == 1  , ...
    der(x) == v , ...
    der(v) == a , ...
    a == -24*(t-0.5)^2, ...
    l == 0 ...
    );
sim = ivp.solve();

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(v) == a, ...
    der(x) == v, ...
    v(t0) == -v(tf) == 1, ...
    x(t0) ==  x(tf) == 0, ...
    x <= l == 1/9 ...
    );
sol = ocp.solve('intervals', 40, 'guess', sim);

figure(1);
subplot(311); hold on
sim.plot(t, x);
sol.plot(t, x);
subplot(312); hold on
sim.plot(t, v);
sol.plot(t, v);
subplot(313); hold on
sim.plot(t, a);
sol.stairs(t, a);


















