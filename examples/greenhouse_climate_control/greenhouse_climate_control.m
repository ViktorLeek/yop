[t, t0, tf] = yop.time();
x = yop.state('x', 2);
u = yop.control();

T = 48;
p = [7.5e-8; 1; 0.1; 4.55e-4; 136.4];
I = max(0, 800*sin(4*pi*t/T - 0.65*pi));
T0 = 15 + 10*sin(4*pi*t/T - 0.65*pi);

dx1 = p(1) * I * x(2);
dx2 = p(2) * (T0-x(2)) + p(3)*u;

ocp = yop.ocp('Greenhouse Climate Control');
ocp.min( at(t==tf, -p(5)*x(1)) + int(p(4)*u) );
ocp.st(t0==0, tf==T, der(x)==[dx1;dx2], x(t0)==[0;10], 0<u<10);

[sol, dms] = ocp.present.solve(100, 4);

time = casadi.Function('x', {dms.w}, {vertcat(dms.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {dms.w}, {horzcat(dms.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {dms.w}, {horzcat(dms.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(311); hold on;
plot(t_sol, 1200*x_sol(:,1))
subplot(312); hold on;
plot(t_sol, x_sol(:,2))
subplot(313); hold on;
stairs(t_sol, [u_sol; nan])