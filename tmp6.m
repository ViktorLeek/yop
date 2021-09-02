%% Problem formulation
clc
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
grp = yop.ocp('Goddard''s Rocket Problem');
grp.max( h(tf) );
grp.st( ...
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
grp.build.present();

%% Transcription

nlp = yop.dms(grp, 10, 4);

%%
w0 = ones(size(nlp.w));

prob = struct('f', nlp.J, 'x', nlp.w, 'g', nlp.g);
solver = casadi.nlpsol('solver', 'ipopt', prob);
sol = solver( ...
    'x0', w0, ...
    'lbx', nlp.w_lb, ...
    'ubx', nlp.w_ub, ...
    'lbg', zeros(size(nlp.g)), ...
    'ubg', zeros(size(nlp.g)) ...
    );

%%
time = casadi.Function('x', {nlp.w}, {vertcat(nlp.t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {nlp.w}, {horzcat(nlp.x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {nlp.w}, {horzcat(nlp.u{:})});
u_sol = full(control(sol.x))';

figure(1)
subplot(411); hold on;
plot(t_sol, x_sol(:,1))
subplot(412); hold on;
plot(t_sol, x_sol(:,2))
subplot(413); hold on;
plot(t_sol, x_sol(:,3))
subplot(414); hold on;
stairs(t_sol, [u_sol; nan])






























