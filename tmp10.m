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
    ... Rocket model
    der(h) == v              , ...
    der(v) == (Wf*c-F_D)/m-g , ...
    der(m) == -Wf            , ...
    ... Initial conditions
    h(t0)==0, ...
    v(t0)==0, ...
    m(t0)==m_max, ...
    ... Box constraints
    h >= 0, ...
    v >=0, ...
    m_min <= m  <= m_max, ...
    Wfmin <= Wf <= Wfmax ...
    );

ocp.build().present();

%%
N = 50;
D = 4;
obj = yop.dc(ocp, N, 'legendre', D);
obj.build()

w0 = ones(size(obj.w));
d = obj.diffcon;
nlp.f = -obj.x(end).y(1);
nlp.x = obj.w;
nlp.g = d;
solver = casadi.nlpsol('solver', 'ipopt', nlp);
sol = solver( ...
    'x0', w0, ...
    'lbx', obj.w_lb, ...
    'ubx', obj.w_ub, ...
    'ubg', zeros(size(d)), ...
    'lbg', zeros(size(d)) ...
    );

%%

tt = [];
for n=1:N
    for r=1:D+1
        tt = [tt(:); obj.t{n,r}];
    end
end
tt = [tt(:); obj.t{N+1,1}];
time = casadi.Function('x', {obj.w}, {tt});
tx_sol = full(time(sol.x));

tt = [];
for n=1:N+1
    tt = [tt(:); obj.t{n,1}];
end
time = casadi.Function('x', {obj.w}, {tt});
tu_sol = full(time(sol.x));

xx = [];
for xk=obj.x
    xx = [xx; xk.y'];
end
state = casadi.Function('x', {obj.w}, {xx});
x_sol = full(state(sol.x));

control = casadi.Function('x', {obj.w}, {horzcat(obj.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(411); hold on;
plot(tx_sol, x_sol(:,1))
subplot(412); hold on;
plot(tx_sol, x_sol(:,2))
subplot(413); hold on;
plot(tx_sol, x_sol(:,3))
subplot(414); hold on;
stairs(tu_sol, [u_sol; nan])