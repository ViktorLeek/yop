[t, t0, tf] = yop.time('t');
p = yop.state('p');   % position
s = yop.state('s');   % speed
a = yop.control('a'); % acceleration
l = 1/9;

x = yop.state('x', 2);
p = x(1);
s = x(2);

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 0.5*int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == trolley(x, a), ...
    ...der(s) == a, ... 
    ...der(p) == s, ...
    p(t0) == 0, s(t0) ==  1, ...
    p(tf) == 0, s(tf) == -1, ...
    p <= l ...
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