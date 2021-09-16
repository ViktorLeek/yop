[t,t0,tf] = yop.time();
x = yop.state();
u = yop.control();

ocp = yop.ocp('Isopermetric Constraint');
ocp.min( int(x) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x) == -sin(x) + u, ...
    x(t0)==1, x(tf)==0, ...
    int(u^2) == 10, ...
    -10 <= x <= 10, ...
    -4 <= u <= 4 ...
    );
[sol, dms] = ocp.present.solve(100, 4);

time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(211); hold on;
plot(t_sol, x_sol(:,1))
subplot(212); hold on;
stairs(t_sol, [u_sol; nan])

%%
[t,t0,tf] = yop.time();
x = yop.state('x', 2);
u = yop.control();

ocp = yop.ocp('Isopermetric Constraint');
ocp.min( int(x(1)) );
ocp.st( ...
    t0==0, tf==1, ...
    der(x(1)) == -sin(x(1)) + u, ...
    der(x(2)) == u^2,  ...
    x(t0)==[1;0], ...
    x(tf)==[0;10], ...
    -10 <= x(1) <= 10, ...
    -4 <= u <= 4 ...
    );
ocp.build.present()

[sol, dms] = ocp.present.solve(100, 4);

time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1) % overlay
subplot(211); hold on;
plot(t_sol, x_sol(:,1))
subplot(212); hold on;
stairs(t_sol, [u_sol; nan])

figure(2)
subplot(311); hold on;
plot(t_sol, x_sol(:,1))
subplot(312); hold on;
plot(t_sol, x_sol(:,2))
subplot(313); hold on;
stairs(t_sol, [u_sol; nan])