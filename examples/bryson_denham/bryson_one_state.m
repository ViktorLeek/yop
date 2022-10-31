yops Times: t t0 tf 
yops State: x size: [2,1]
yops Controls: u

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(u^2) );
ocp.st( 0 == t0 < tf == 1 );
ocp.st( der(x) == [x(2); u] );
ocp.st(  x(t0) == [0; +1] );
ocp.st(  x(tf) == [0; -1] );
ocp.st(  x(1)  <= 1/9     );
sol = ocp.solve('ival', 20);

figure(1);
subplot(211); hold on
sol.plot(t, x(2));
subplot(212); hold on
sol.stairs(t, u);

%% Bryson-Denham, 1 state

% Approximate abs with a smooth approximation
k = 300; % 700 - limit
abs = @(x) 2/k*log(1+exp(k*x)) - x - 2/k*log(2);

t = yop.time();
t0 = yop.time0();
tf = yop.timef();
u = yop.control();
v = yop.state();

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(der(v)^2) );
ocp.st( t0 == 0, tf == 1 );
ocp.st( der(v) == u );
ocp.st(  v(t0) == +1 );
ocp.st(  v(tf) == -1 );
ocp.st( int(v) == 0 ); % x(0) == x(1) == 0
ocp.st( int( abs(v) ) <= 2*1/9 ); % x <= 1/9
sol = ocp.solve('ival', 20);

figure(1);
subplot(211); hold on
sol.plot(t, v);
subplot(212); hold on
sol.stairs(t, u);
