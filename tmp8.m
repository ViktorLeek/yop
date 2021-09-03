[t, t0, tf] = yop.time('t');
p = yop.state('p'); % position
s = yop.state('s'); % speed
a = yop.control('a'); % acceleration
l = 1/9;

ocp = yop.ocp('Bryson-Denham Problem');
ocp.min ( 0.5*int(a^2) );
ocp.st( ...
    t0==0, tf==1, ...
    der(p) == s, ...
    der(s) == a, ...
    der(i) == a^2, ...
    p(t0) == 0, s(t0) ==  1, ...
    p(tf) == 0, s(tf) == -1, ...
    p <= l ...
    );
ocp.build().present();

%%
nlp = yop.dms(ocp, 20, 10);

w0 = ones(size(nlp.w));

g = vertcat(nlp.eq{:});
h = vertcat(nlp.ieq{:});
prob = struct('f', nlp.J, 'x', nlp.w, 'g', [g; h]);
solver = casadi.nlpsol('solver', 'ipopt', prob);
sol = solver( ...
    'x0', w0, ...
    'lbx', nlp.w_lb, ...
    'ubx', nlp.w_ub, ...
    'ubg', [zeros(size(g)); zeros(size(h))], ...
    'lbg', [zeros(size(g)); -inf(size(h))] ...
    );

%%
time = casadi.Function('x', {nlp.w}, {vertcat(nlp.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {nlp.w}, {horzcat(nlp.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {nlp.w}, {horzcat(nlp.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(311); hold on;
plot(t_sol, x_sol(:,1))
subplot(312); hold on;
plot(t_sol, x_sol(:,2))
subplot(313); hold on;
stairs(t_sol, [u_sol; nan])