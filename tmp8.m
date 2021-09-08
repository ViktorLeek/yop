[t, t0, tf] = yop.time('t');
x = yop.state('x'); % position
v = yop.state('v'); % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( int(0.5*int(a^2) + 0*x(t==0.5)) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == v, ...
    der(v) == a, ...
    x(t0) == 0, v(t0) ==  1, ...
    x(tf) == 0, v(tf) == -1, ...
    x <= l ...
    ...x(t==0.5) <= 0.1 ...
    ... x(t==0.3) <= 0.09 ...
    );
ocp.build().present();
%%
dms = yop.dms(ocp, 40, 4);
sol = dms.solve();

%%
time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(311); hold on;
plot(t_sol, x_sol(:,1))
subplot(312); hold on;
plot(t_sol, x_sol(:,2))
subplot(313); hold on;
stairs(t_sol, [u_sol; nan])