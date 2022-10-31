t  = yop.time();
t0 = yop.time0();
tf = yop.timef();
s  = yop.state();
v  = yop.state();
u  = yop.control();

s = @(tp) s(t==tp);
v = @(tp) v(t==tp);
u = @(tp) u(t==tp);

ocp = yop.ocp();
ocp.min(1/2 * int(u(t)^2));
ocp.st(t0==0, tf==1, ...
       der(s(t)) == v(t), ...
       der(v(t)) == u(t), ...
       s(0) ==  s(1) == 0, ...
       v(0) == -v(1) == 1, ...
       s(t) <= 1/9 );

sol = ocp.solve();
figure(1); hold on;
sol.plot(t, s(t));
sol.plot(t, v(t));
sol.stairs(t, u(t)/5);