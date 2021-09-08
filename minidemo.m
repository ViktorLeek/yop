%% Variations on Bryson-Denham problem
% 1) Basic problem
[t, t0, tf] = yop.time('t');
x = yop.state('x'); % position
v = yop.state('v'); % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( 0.5*int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == v, ...
    der(v) == a, ...
    x(t0) == 0, ...
    v(t0) ==  1, ...
    x(tf) == 0, ...
    v(tf) == -1, ...
    x <= l ...
    );
ocp.build().present();

dms = yop.dms(ocp, 40, 4);
sol = dms.solve();
plot_res(dms, sol);

%% Using function to hold differential eq. - better for large systems
clear
[t, t0, tf] = yop.time('t');
x = yop.state('x', 2); % [position; speed]
u = yop.control('u');
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( 0.5*int(u^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == trolley(x,u), ...
    x(t0) == [0; 1], ...
    x(tf) ==  [0; -1], ...
    x(1) <= l ...
    );
ocp.build().present();

dms = yop.dms(ocp, 40, 4);
sol = dms.solve();
plot_res(dms, sol);

%% Scalar states with vector differential eq. - Good for large systems with
%  complex constraints, as 
clear
[t, t0, tf] = yop.time('t');
pos = yop.state('x'); % position
speed = yop.state('v'); % speed
acc = yop.control('a'); % acceleration
l = 1/9;

% Vectorize
x_vec = [pos; speed];
u_vec = acc;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( 0.5*int(acc^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x_vec) == trolley(x_vec, u_vec), ...
    pos(t0) == 0, ...
    speed(t0) ==  1, ...
    pos(tf) == 0, ...
    speed(tf) == -1, ...
    pos <= l ...
    );
ocp.build().present();

dms = yop.dms(ocp, 40, 4);
sol = dms.solve();
plot_res(dms, sol);

%% First formulation with complex constraints
clear
% 1) Basic problem
[t, t0, tf] = yop.time('t');
x = yop.state('x'); % position
v = yop.state('v'); % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( 0.5*int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == v, ...
    der(v) == a, ...
    x(t0) == 0, ...
    v(t0) ==  1, ...
    x(tf) == 0, ...
    v(tf) == -1, ...
    x <= l, ...
    x(t==0.5) <= 0.1 ... Evaluate x at arbitrary timepoint
    );
ocp.build().present();

dms = yop.dms(ocp, 40, 4);
sol = dms.solve();
plot_res(dms, sol);

%% First formulation with complex objective
% 1) Basic problem
[t, t0, tf] = yop.time('t');
x = yop.state('x');   % position
v = yop.state('v');   % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min( int( 0.5*int(a^2) - 10*v(t==0.5)*v^2 ) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == v, ...
    der(v) == a, ...
    x(t0) == 0, ...
    v(t0) ==  1, ...
    x(tf) == 0, ...
    v(tf) == -1, ...
    ... x(t==0.5) <= 0.1, ... Evaluate x at arbitrary timepoint
    x <= l ...
    );
ocp.build().present();

dms = yop.dms(ocp, 40, 4);
sol = dms.solve();
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