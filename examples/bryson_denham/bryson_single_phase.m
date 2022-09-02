%% Original formulation
yops Times: t t0 tf % Parsed by position: t, t0, tf
yops States: x v    % position, speed
yops Control: a     % acceleration
yops Param: l       % maximum cart position

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
sol.stairs(t, a);

%% Clean implementation
yops Times: t t0 tf 
yops States: x size: [2,1]
yops Controls: u

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(u^2) );
ocp.st( 0 == t0 < tf == 1 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve();

figure(1);
subplot(211); hold on
sol.plot(t, x);
subplot(212); hold on
sol.stairs(t, u);

%% State vector, minimum value and traveled distance
%   Piecewise quadratic control input (deg == 2)
yops Times: t t0 tf 
yops States: x size: [2,1]
yops Controls: u int: 2 % picewise quadratic control input

J = 1/2 * int(u^2);
ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( J );
ocp.st(   t0==0, tf==1    );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve('ival', 50);

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
yops Times: t t0 tf States: x size: [2,1] Controls: u

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(u^2) );
ocp.st( t0==0, tf==1, ...
    der(x) == [x(2); u], ...
     x(t0) == [0; 1], ...
     x(tf) == [0; -1], ...
     x(1)  <= 1/9 ...
    );

sol = ocp.solve('ival', 15, 'dx', 2);
figure(1);
subplot(311); hold on
sol.plot(t, x(1), 'mag', 5);
subplot(312); hold on
sol.plot(t, x(2));
subplot(313); hold on
sol.stairs(t, u);

%% Compact representation
[t0, tf, t, x, u] = yop.vars('nx', 2, 'nu', 1);
ocp = yop.ocp().min( 1/2 * int(u^2) );
ocp.st(   t0==0, tf==1    );
ocp.st( der(x) == [x(2); u]  );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve('ival', 20, 'dx', 2);
figure(1); hold on
sol.plot(t, x);
sol.stairs(t, u);

%% Minreal
[t0, tf, t, x, u] = yop.vars('nx', 2, 'nu', 1);
yop.ocp().min(1/2*int(u^2)).st(t0==0,tf==1,der(x)==[x(2);u],x(t0)==[0;1],...
    x(tf)==[0;-1], x(1)<=1/9).solve().plot(t, [x;u]);

%% Trade-off between control effort and traveled distance
yops Times: t t0 tf States: x v Control: a Param: l

beta = 30; % trade-off factor

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) + beta*l );
ocp.st( t0==0, tf==1 );
ocp.st( der(v) == a );
ocp.st( der(x) == v );
ocp.st( v(t0) == -v(tf) == 1 );
ocp.st( x(t0) == x(tf) == 0 );
ocp.st( x <= l );

sol = ocp.solve('ival', 20, 'dx', 2);

figure(1);
subplot(311); hold on
sol.plot(t, x, 'mag', 5);
sol.plot(t, 0*t + l);
subplot(312); hold on
sol.plot(t, v, 'mag', 5);
subplot(313); hold on
sol.stairs(t, a, 'mag', 5);

%% Simulated initial guess
yops Times: t t0 tf 
yops States: x v 
yops Control: a 
yops Parameter: l

sim = yop.ivp( ...
    t0==0, tf==1, ...
    x(t0) == 0  , ...
    v(t0) == 1  , ...
    der(x) == v , ...
    der(v) == a , ...
    a == -24*(t-0.5)^2, ...
    l == 0 ...
    );
res = sim.solve();

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
sol = ocp.solve('ival', 40, 'guess', res);

figure(1);
subplot(311); hold on
res.plot(t, x);
sol.plot(t, x);
subplot(312); hold on
res.plot(t, v);
sol.plot(t, v);
subplot(313); hold on
res.plot(t, a);
sol.stairs(t, a);

%% Time trade-off
yops Times: t t0 tf % Parsed by position: t, t0, tf
yops State: x v     % position, speed
yops Ctrls: a       % acceleration
yops Param: l       % maximum cart position

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) + 0.1*tf );
ocp.st( t0==0, tf==1 );
ocp.st( der(v) == a );
ocp.st( der(x) == v );
ocp.st( v(t0) == -v(tf) == 1 );
ocp.st( x(t0) ==  x(tf) == 0 );
ocp.st( x <= l == 1/9 );
sol = ocp.solve('ival', 10, 'dx', 2);

figure(1);
subplot(311); hold on
sol.plot(t, x);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.stairs(t, a);

%% Using (t==0) and (t==1) instead of (t0) and (tf)
t  = yop.independent();
t0 = yop.time0();
tf = yop.timef();
s  = yop.state();
v  = yop.state();
u  = yop.control();

ocp = yop.ocp();
ocp.min( 1/2 * int(u^2) );
ocp.st( t0==0, tf==1 );
ocp.st( der(s) == v );
ocp.st( der(v) == u );
ocp.st( s(t==0) ==  s(t==1) == 0 );
ocp.st( v(t==0) == -v(t==1) == 1 );
ocp.st( s <= 1/9 );

sol = ocp.solve();
figure(1); hold on;
sol.plot(t, s);
sol.plot(t, v);
sol.plot(t, u);
















