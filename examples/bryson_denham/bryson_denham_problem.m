%% Original formulation
[t, t0, tf] = yop.time('t');
p = yop.state('p'); % position
s = yop.state('s'); % speed
a = yop.control('a'); % acceleration
l = yop.parameter('l'); % maximum position of the cart

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(s) == a, ...
    der(p) == s, ...
    p(t0) == p(tf) == 0, ...
    s(t0) == -s(tf) == 1, ...
    p <= l == 1/9 ... An interesting variation is to balance l and control effort
    );

[tt,xx,uu,pp]=ocp.solve('method', 'dms', 'intervals', 20, 'rk4_steps', 4);

figure(1)
subplot(311); hold on;
plot(tt, xx(:,1))
subplot(312); hold on;
plot(tt, xx(:,2))
subplot(313); hold on;
stairs(tt, [uu; nan])
%% Guaranteed box constraints for boundary conditions

[t, t0, tf] = yop.time('t');
p = yop.state('p'); % position
s = yop.state('s'); % speed
a = yop.control('a'); % acceleration
l = yop.parameter('l'); % maximum position of the cart

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(p) == s, ...
    der(s) == a, ...
    p(t0) == 0 == p(tf), ... % SRF -> {p(t0)==1, 1==p(tf)}
    s(t0) ==  1, ...
    s(tf) == -1, ...
    p <= l == 1/9 ...
    );

[tt,xx,uu,pp]=ocp.solve('method', 'dms', 'intervals', 20, 'rk4_steps', 4);

figure(1)
subplot(311); hold on;
plot(tt, xx(:,1))
subplot(312); hold on;
plot(tt, xx(:,2))
subplot(313); hold on;
stairs(tt, [uu; nan])


%% Removal of unncessary parameter

[t, t0, tf] = yop.time('t');
p = yop.state('p'); % position
s = yop.state('s'); % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(p) == s, ...
    der(s) == a, ...
    p(t0) == 0, s(t0) ==  1, ...
    p(tf) == 0, s(tf) == -1, ...
    p <= l ...
    );
ocp.present();
[sol, dms] = ocp.present.solve(20, 4);
plot_res(dms, sol);


%%

[t, t0, tf] = yop.time('t');
x = yop.state('x', 2); % position
u = yop.control('u'); % acceleration

p = x(1);
s = x(2);
a = u;

l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( 1/2 * int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == [s; a], ...
    p(t0) == 0, s(t0) ==  1, ...
    p(tf) == 0, s(tf) == -1, ...
    p(0.2 < t < 0.8) <= 1/10, ...
    p <= l ...
    );
ocp.present();
[sol, dms] = ocp.present.solve(20, 4);
plot_res(dms, sol);

%%
function plot_res(dms, sol)
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
end