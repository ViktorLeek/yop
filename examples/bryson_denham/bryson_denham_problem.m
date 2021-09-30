%% Original formulation
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');
x  = yop.state('x');     % position
v  = yop.state('v');     % speed
a  = yop.control('a');   % acceleration
l  = yop.parameter('l'); % maximum position of the cart

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) + 30*l );
ocp.st( ...
    t0==0, tf==1, ...
    der(v) == a, ...
    der(x) == v, ...
    v(t0) == -v(tf) == 1, ...
    x(t0) == x(tf) == 0, ...
    x <= l == 1/9 ...
    );

sol = ocp.solve('intervals', 20);

figure(1);
subplot(311); hold on
sol.plot(t, x);
subplot(312); hold on
sol.plot(t, v);
subplot(313); hold on
sol.stairs(t, a);

%% Guaranteed box constraints for boundary conditions
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');
x  = yop.state('x');     % position
v  = yop.state('v');     % speed
a  = yop.control('a');   % acceleration

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
    p == int(a^2), ...
    x <= 1/9 ...
    );

sol = ocp.solve('intervals', 20);

figure(1);
subplot(311); hold on
sol.plot(t, x, 'refine', 5);
subplot(312); hold on
sol.plot(t, v, 'refine', 5);
subplot(313); hold on
sol.stairs(t, a, 'refine', 5);


%% Compact representation
[t0, tf, t, x, u] = yop.ocp_variables('nx', 2, 'nu', 1);
ocp = yop.ocp().min( 1/2 * int(u^2) );
ocp.st(tf==1, ...
    der(x) == [x(2); u], ...
    x(t0)  == [0;  1], ...
    x(tf)  == [0; -1], ...
    x(1)   <= 1/9 ...
    );
sol = ocp.solve('intervals', 20);
figure(1);
subplot(211); hold on
sol.plot(t, x);
subplot(212); hold on
sol.stairs(t, u);

%% Minreal
[t0, tf, t, x, u] = yop.ocp_variables('nx', 2, 'nu', 1);
yop.ocp().min(1/2*int(u^2)).st(tf==1, der(x)==[x(2);u], x(t0)==[0; 1], ...
    x(tf)==[0;-1], x(1)<=1/9).solve('intervals',20).plot(t,[x;u]);

%% Trade-off between control effort and traveled distance
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');
x  = yop.state('x'); % position
v  = yop.state('v'); % speed
a  = yop.control('a'); % acceleration
l  = yop.parameter('l'); % maximum position of the cart

beta = 30; % trade-off factor

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) + beta*l );
ocp.st( ...
    t0==0, tf==1, ...
    der(v) == a, ...
    der(x) == v, ...
    v(t0) == -v(tf) == 1, ...
    x(t0) == x(tf) == 0, ...
    x <= l ...
    );

sol = ocp.solve('intervals', 20);

figure(1);
subplot(311); hold on
sol.plot(t, x, 'refine', 5);
sol.plot(t, 0*t + l);
subplot(312); hold on
sol.plot(t, v, 'refine', 5);
subplot(313); hold on
sol.stairs(t, a, 'refine', 5);
