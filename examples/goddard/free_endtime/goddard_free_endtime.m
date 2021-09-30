%% Formulation 1
[t, t0, tf] = yop.time('t');
x = yop.state('x', 3);
u = yop.control('u');

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( rocket.height(tf) );
% ocp.max( int(rocket.velocity) );
ocp.st( ...
    t > 0, ...
    t0==0, ...tf==212, ...
    der(x) == rocket_model(x, u), ...
    rocket.height(t0)   == 0    , ...
    rocket.velocity(t0) == 0    , ...
    rocket.mass(t0)     == 215  , ...
    68 <= rocket.mass <= 215, ...
    0 <= rocket.fuel_mass_flow <= 9.5 ...
    );

sol = ocp.solve('intervals', 50);

%%
figure(1);
subplot(411); hold on
sol.plot(t, x(1));
sol.plot(t(t==1), x(1).at(t==1), 'x')
plot(sol.resolved_value(t, 1000), sol.resolved_value(x(1), 1000))
subplot(412); hold on
sol.plot(t, x(2));
subplot(413); hold on
sol.plot(t, x(3));
subplot(414); hold on
sol.stairs(t, u);

%%
expr = yop.ocp_expr(int(x(1)));

[vars, tps, ints, ders, sn] = yop.ocp.find_special_nodes(expr.ast);

args = { ...
    mx_vec(ocp.independent), ...
    mx_vec(ocp.states), ...
    mx_vec(ocp.controls), ...
    mx_vec(ocp.parameters), ...
    mx_vec(tps), ...
    mx_vec(ints), ...
    mx_vec(ders) ...
    };

set_mx(ocp.variables);
set_mx([tps, ints, ders]);

for node = [tps, ints, ders]
    mx_expr = fw_eval(node.ast.expr);
    node.fn = casadi.Function('fn', args, {mx_expr});
end

expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});

t0 = full(sol.x(1));
tf = full(sol.x(2));
[tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), n_elem(ints), ...
    nlp.N, nlp.tau, nlp.dt, t0, tf, nlp.t, nlp.x, nlp.u, nlp.p);

disc = yop.parameterize_expression( ...
    expr, ...
    nlp.N, ...
    nlp.tau, ...
    nlp.t, ...
    nlp.x, ...
    nlp.u, ...
    nlp.p, ...
    tpv, ...
    intv, ...
    [] ...
    );

f = casadi.Function('f', {nlp.w}, {disc});
num = full(f(sol.x))

%%

time = casadi.Function('t', {nlp.w}, {mat(nlp.t)});
tt = full(time(sol.x));

state = casadi.Function('x', {nlp.w}, {mat(nlp.x)});
xx = full(state(sol.x));

control = casadi.Function('u', {nlp.w}, {mat(nlp.u)});
uu = full(control(sol.x));

parameter = casadi.Function('p', {nlp.w}, {nlp.p});
pp = full(parameter(sol.x));

tlp = yop.collocated_time(tt(1), tt(end), nlp.N);

y = cell(nlp.N+1,1);
for n=1:nlp.N
    y{n} = xx(:, (nlp.d+1)*n-nlp.d:(nlp.d+1)*n);
end
y{nlp.N+1} = xx(:,end);
xlp = yop.collocated_expression(nlp.N, nlp.tau, y);

y = cell(nlp.N+1,1);
for n=1:nlp.N
    y{n} = uu(:, n);
end
y{nlp.N+1} = uu(:,end);
ulp = yop.collocated_expression(nlp.N, 0, y);

% expr = yop.ocp_expr([t; x; u]);
expr = yop.ocp_expr(x(1).at(t==212));

[vars, tps, ints, ders, sn] = yop.ocp.find_special_nodes(expr.ast);

args = { ...
    mx_vec(ocp.independent), ...
    mx_vec(ocp.states), ...
    mx_vec(ocp.controls), ...
    mx_vec(ocp.parameters), ...
    mx_vec(tps), ...
    mx_vec(ints), ...
    mx_vec(ders) ...
    };

set_mx(ocp.variables);
set_mx([tps, ints, ders]);

for node = [tps, ints, ders]
    mx_expr = fw_eval(node.ast.expr);
    node.fn = casadi.Function('fn', args, {mx_expr});
end

expr.fn = casadi.Function('fn', args, {fw_eval(expr.ast)});

[tpv, intv] = yop.param_special_nodes(sn, n_elem(tps), n_elem(ints), ...
    nlp.N, nlp.tau, nlp.dt, full(sol.x(1)), full(sol.x(2)), nlp.t, nlp.x, nlp.u, nlp.p);

tpf = casadi.Function('f', {nlp.w}, {tpv});
intf = casadi.Function('f', {nlp.w}, {intv});

TP = tpf(sol.x);
I = intf(sol.x);

if expr.is_transcription_invariant
    v = expr.fn( ...
        tlp(1).evaluate(0), ...
        xlp(1).evaluate(0), ...
        ulp(1).evaluate(0), ...
        pp, TP, I, []);
else
    pnts = 200;
    grid = linspace(tlp(1).evaluate(0), tlp(end).evaluate(0), pnts);
    v = [];
    for tp=grid
        vk = expr.fn( ...
            tlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
            xlp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
            ulp.value(tp, tlp(1).evaluate(0), tlp(end).evaluate(0)), ...
            pp, TP, I, []);
        v = [v, vk];
    end
end
V = full(v);
%%
figure(1);
subplot(411); hold on
plot(V(1,:), V(2,:))
subplot(412); hold on
plot(V(1,:), V(3,:))
subplot(413); hold on
plot(V(1,:), V(4,:))
subplot(414); hold on
stairs(V(1,:), V(5,:))

%%
figure(1);
subplot(411); hold on
plot(V(1,:), V(2,:), 'x')
subplot(412); hold on
plot(V(1,:), V(3,:), 'x')
subplot(413); hold on
plot(V(1,:), V(4,:), 'x')
subplot(414); hold on
stairs(V(1,:), V(5,:), 'x')

% plot(sol.value(t), sol.value(x))
% subplot(212); hold on
% stairs(sol.value(t), sol.value(u))
%%
[tt, xx, uu, pp, tx] = ocp.solve('method', 'dc', 'intervals', 50);

figure(1)
subplot(411); hold on;
plot(tx, xx(:,1))
subplot(412); hold on;
plot(tx, xx(:,2))
subplot(413); hold on;
plot(tx, xx(:,3))
subplot(414); hold on;
stairs(tt, [uu; nan])

%% Formulation 1 DMS
[t, t0, tf] = yop.time('t');
x = yop.state('x', 3);
u = yop.control('u');

[~, y] = rocket_model(x, u);
rocket = y.rocket;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( int(rocket.velocity) );
ocp.st( ...
    t0==0, ...
    ... dynamics 
    p == int(rocket.velocity), ...
    der(x) == rocket_model(x, u), ...
    ... initial value
    rocket.height(t0)   == 0, ...
    rocket.velocity(t0) == 0, ...
    rocket.mass(t0)     == 215, ...
    ... Box constraints
    68 <= rocket.mass <= 215, ...
     0 <= rocket.fuel_mass_flow <= 9.5 ...
    );
[tt, xx, uu, pp] = ocp.solve('method', 'dms', 'intervals', 50);

figure(1)
subplot(411); hold on;
plot(tt, xx(:,1))
subplot(412); hold on;
plot(tt, xx(:,2))
subplot(413); hold on;
plot(tt, xx(:,3))
subplot(414); hold on;
stairs(tt, [uu; nan])

%% Formulation 2
[t, t0, tf] = yop.time('t');
v = yop.state('v');
h = yop.state('h');
m = yop.state('m');
x = [v; h; m];

Wf = yop.control('Wf');
u = Wf;

m_max = 215;
m_min = 68;
Wf_max = 9.5;

ocp = yop.ocp('Goddard''s Rocket Problem');
ocp.max( h(tf) );
ocp.st( ...
    t0 == 0, ...
    der(x) == rocket_model(x, u), ...    
    v(t0) == 0, ...
    h(t0) == 0, ...
    m(t0) == m_max, ...
    h >= 0, ...
    v >= 0, ...
    m_min <= m <= m_max, ...
    0 <= Wf <= Wf_max ...
    );
ocp.build();


[tt, xx, uu, pp] = ocp.solve('method', 'dms', 'intervals', 70);

figure(1)
subplot(411); hold on;
plot(tt, xx(:,1))
subplot(412); hold on;
plot(tt, xx(:,2))
subplot(413); hold on;
plot(tt, xx(:,3))
subplot(414); hold on;
stairs(tt, [uu; nan])

%% Formulation 3
% Time
[t, t0, tf] = yop.time('t');

% Rocket parameters
D0 = 0.01227; beta = 0.145e-3; c = 2060;

% Constants
g0 = 9.81; r0 = 6.371e6;

% Rocket model
h  = yop.state('h');  % Rocket height
v  = yop.state('v');  % Rocket speed
m  = yop.state('m');  % Rocket mass
Wf = yop.control('Wf');  % Rocket fuel massflow

% Drag force and gravitational acceleration
F_D = D0 * exp(-beta*h) * v^2;
g   = g0*(r0/(r0+h))^2;

% Mass boundaries
m_min = 68; m_max = 215; 

% Control boundaries
Wfmin = 0; Wfmax = 9.5;

% Optimal control problem
ocp = yop.ocp();
ocp.max( h(tf) );
ocp.st( ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_max, ...
    ... Rocket model
    der(v) == (Wf*c-F_D)/m-g , ...
    der(h) == v              , ...
    der(m) == -Wf            , ...
    ... Box constraints
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

[tt, xx, uu, pp, tx] = ocp.solve('method', 'dc', 'intervals', 50);

figure(1)
subplot(412); hold on;
plot(tx, xx(:,1))
subplot(411); hold on;
plot(tx, xx(:,2))
subplot(413); hold on;
plot(tx, xx(:,3))
subplot(414); hold on;
stairs(tt, [uu; nan])


















