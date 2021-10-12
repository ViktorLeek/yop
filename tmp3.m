%% Original formulation
t0 = yop.time0('t0');
tf = yop.timef('tf');
t  = yop.time('t');
x  = yop.state('x');     % position
v  = yop.state('v');     % speed
a  = yop.control('a');   % acceleration
l  = yop.parameter('l'); % maximum position of the cart

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
sol = ocp.solve();























