%%
% syms tt xx1 xx2 xx3 uu
% xx = [xx1; xx2; xx3];

import yop.*
[t, t0, tf] = independent('t');
x = state('x', 3);
u = control('u');

[dx, y] = rocket_model(x, u);

m0 = 215;
mf = 68;

ocp = yop.ocp();
ocp.max(y.rocket.height(tf));
ocp.st(...
    ... dynamics 
    der(x) == rocket_model(x, u), ...
    ... initial condition
    t0  == 0, ...
    y.rocket.height(t0) == 0, ...
    y.rocket.velocity(t0) == 0, ...
    y.rocket.mass(t0) == m0, ...
    ... Box constraints
    y.rocket.height >= 0, ...
    y.rocket.velocity >= 0, ...
    mf <= y.rocket.mass <= m0, ...
     0 <= y.rocket.fuel_mass_flow <= 9.5 ...
    );
ocp.build();

% x.m_value = xx;
% u.m_value = uu;

% rocket.height(t0);
% rocket.height(tf);
% rocket.height(t==t0);
% rocket.height(t==tf);
% rocket.height(t==4);

%% Multiple shooting

% Separate the problem variables into categories
var = struct;
var.t = {};
var.t0 = {};
var.tf = {};
var.x = {};
var.z = {};
var.u = {};
var.p = {};
for k=1:length(ocp.variables)
    switch class(ocp.variables{k})
        case 'yop.ast_independent'
            var.t = {var.t{:}, ocp.variables{k}};
            
        case 'yop.ast_independent_initial'
            var.t0 = {var.t0{:}, ocp.variables{k}};
            
        case 'yop.ast_independent_final'
            var.tf = {var.tf{:}, ocp.variables{k}};
            
        case 'yop.ast_state'
            var.x = {var.x{:}, ocp.variables{k}};
            
        case 'yop.ast_algebraic'
            var.z = {var.z{:}, ocp.variables{k}};
            
        case 'yop.ast_control'
            var.u = {var.u{:}, ocp.variables{k}};
            
        case 'yop.ast_parameter'
            var.p = {var.p{:}, ocp.variables{k}};
            
        otherwise
            error('[yop] Error: Unknown variable type.')
    end
end

% Populate the variables values
sym = struct;
sym.t = {};
sym.t0 = {};
sym.tf = {};
sym.x = {};
sym.z = {};
sym.u = {};
sym.p = {};

sym.t = casadi.MX.sym('t');
for k=1:length(var.t)
    var.t{k}.m_value = sym.t;
end

sym.t0 = casadi.MX.sym('t0');
for k=1:length(var.t0)
    var.t0{k}.m_value = sym.t0;
end

sym.tf = casadi.MX.sym('tf');
for k=1:length(var.tf)
    var.tf{k}.m_value = sym.tf;
end

for k=1:length(var.x)
    sz = size(var.x{k});
    xk = casadi.MX.sym('x', sz(1), sz(2));
    sym.x{k} = xk;
    var.x{k}.m_value = xk;
end

for k=1:length(var.z)
    sz = size(var.z{k});
    zk = casadi.MX.sym('z', sz(1), sz(2));
    sym.z{k} = zk;
    var.z{k}.m_value = zk;
end

for k=1:length(var.u)
    sz = size(var.u{k});
    uk = casadi.MX.sym('u', sz(1), sz(2));
    sym.u{k} = uk;
    var.u{k}.m_value = uk;
end

for k=1:length(var.p)
    sz = size(var.p{k});
    pk = casadi.MX.sym('p', sz(1), sz(2));
    sym.p{k} = pk;
    var.p{k}.m_value = pk;
end

t = sym.t;
t0 = sym.t0;
tf = sym.tf;
x = vertcat(sym.x{:}); % Could lead to problems if dimensions are wrong
z = vertcat(sym.z{:});
u = vertcat(sym.u{:});
p = vertcat(sym.p{:});

inputs = {t, t0, tf, x, z, u, p};
% inputs = inputs(~cellfun('isempty', inputs));

%% Temporary

t  = casadi.MX.sym('t');
t0 = casadi.MX.sym('t0');
tf = casadi.MX.sym('tf');
x  = casadi.MX.sym('x', 3);
z  = casadi.MX.sym('z', 0);
u  = casadi.MX.sym('u');
p  = casadi.MX.sym('p', 0);
inputs = {t, t0, tf, x, u, p};

var.x{1}.m_value = x;
var.u{1}.m_value = u;
var.t0{1}.m_value = t0;

%% Dynamics

% As the variables node the nodes have been initialized with casadi
% variables as values, the expressions are evaluated, and at the root the
% relevant expression is obtained.

% This should be done for every element of ocp.differential
[sort, n_elem] = topological_sort(ocp.differential{1}.expr);

for k=1:(n_elem-1)
    forward(sort{k});
end
xdot_expr = forward(sort{n_elem}); % root holds the dynamics.
f = casadi.Function('xdot', inputs, {xdot_expr, 0}); % zero is for objective

%% Objective function
[sort, n_elem] = topological_sort(ocp.objective);

for k=1:(n_elem-1)
    forward(sort{k});
end
J_expr = forward(sort{n_elem}); % root holds the objective
J_fn = casadi.Function('J', inputs, {J_expr});


%% Integrator
T = tf-t0; % Time horizon
N = 50; % number of control intervals

% Fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = T/N/M;
T0  = casadi.MX.sym('T');
X0 = casadi.MX.sym('X0', size(x,1));
U  = casadi.MX.sym('U', size(u,1));
P  = casadi.MX.sym('P', size(p,1));
T = T0;
X = X0;
Q = 0;
for j=1:M
    [k1, k1_q] = f(T       , t0, tf, X            , U, P);
    [k2, k2_q] = f(T + DT/2, t0, tf, X + DT/2 * k1, U, P);
    [k3, k3_q] = f(T + DT/2, t0, tf, X + DT/2 * k2, U, P);
    [k4, k4_q] = f(T + DT  , t0, tf, X + DT   * k3, U, P);
    T = T + DT;
    X = X + DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = casadi.Function('F', {T0, t0, tf, X0, U, P}, {X, Q});

%% Constraints

t0_lb = 0;
t0_ub = 0;

tf_lb = 0;
tf_ub = inf;

x0_lb = [0; 0; m0];
x0_ub = [0; 0; m0];

x_lb = [0; 0; mf];
x_ub = [inf; inf; m0];

xf_lb = [0; 0; mf];
xf_ub = [inf; inf; m0];

u_lb = 0;
u_ub = 9.5;

%% Initial guess
xg = ones(3,1);
ug = 1;

%% Variables NLP
nx = 3;
nu = 1;
np = 0;

% Discretized independent
t0 = casadi.MX.sym('t0');
tf = casadi.MX.sym('tf');

% Discretized state
x = cell(N+1, 1);
for k=1:(N+1)
    x{k} = casadi.MX.sym(['x_' num2str(k)], nx);
end

% Discretized control
u = cell(N, 1);
for k=1:N
    u{k} = casadi.MX.sym(['u_' num2str(k)], nu);
end

% Parameters
p = casadi.MX.sym('p', np);

% Compute grid timepoints
t = cell(N+1, 1);
dt = (tf-t0)/N;
t{1} = t0;
t{end} = tf;
for k=2:N
    t{k} = t{k-1} + dt;
end

% Integrate objective
J = 0;
for k=1:N
    [~, qf] = F(t{k}, t0, tf, x{k}, u{k}, p);
    J = J + qf;
end
J = J + -x{end}(2);

% Integrate dynamics
xf = cell(N, 1);
for k=1:N
    xf{k} = F(t{k}, t0, tf, x{k}, u{k}, p);
end

% Introduce defect constraints
g_dyn = cell(N, 1);
for k=1:N
    g_dyn{k} = xf{k} - x{k+1};
end
g_dyn_lb = zeros(N*nx, 1);
g_dyn_ub = zeros(N*nx, 1);

% Box constraints
t0_lb = t0_lb;
t0_ub = t0_ub;

tf_lb = tf_lb;
tf_ub = tf_ub;

x_lb = repmat(x_lb, N+1, 1);
x_ub = repmat(x_ub, N+1, 1);
x_lb(1:3) = x0_lb;
x_ub(1:3) = x0_ub;
x_lb(end-2:end) = xf_lb;
x_ub(end-2:end) = xf_ub;

u_lb = repmat(u_lb, N, 1);
u_ub = repmat(u_ub, N, 1);

p_lb = -inf(np, 1);
p_ub =  inf(np, 1);

w = vertcat(t0, tf, x{:}, u{:}, p);
w_lb = [t0_lb; tf_lb; x_lb; u_lb; p_lb];
w_ub = [t0_ub; tf_ub; x_ub; u_ub; p_ub];
w0 = ones(size(w));

prob = struct('f', J, 'x', w, 'g', vertcat(g_dyn{:}));
solver = casadi.nlpsol('solver', 'ipopt', prob);
sol = solver('x0', w0, 'lbx', w_lb, 'ubx', w_ub, 'lbg', g_dyn_lb, 'ubg', g_dyn_ub);


%% Plot solution

time = casadi.Function('x', {w}, {vertcat(t{:})});
t_sol = full(time(sol.x));

state = casadi.Function('x', {w}, {horzcat(x{:})});
x_sol = full(state(sol.x))';

control = casadi.Function('x', {w}, {horzcat(u{:})});
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


























%% Formulate NLP
% nx = size(x,1);
% nu = size(u,1);
% np = size(p,1);
% 
% % Start with an empty NLP
% w={};
% w0 = [];
% lbw = [];
% ubw = [];
% J = 0;
% g={};
% lbg = [];
% ubg = [];
% 
% t0 = casadi.MX.sym('t0');
% w = {w{:}, t0};
% lbw = [lbw; t0_lb];
% ubw = [ubw; t0_ub];
% w0 = [w0; 0];
% 
% tf = casadi.MX.sym('tf');
% w = {w{:}, tf};
% lbw = [lbw; tf_lb];
% ubw = [ubw; tf_ub];
% w0 = [w0; 200];
% 
% p = casadi.MX.sym('p', np);
% w = {w{:}, p};
% % lbw = [lbw; 0];
% % ubw = [ubw; 0];
% % w0 = [w0; 0];
% 
% % "Lift" initial conditions
% Xk = casadi.MX.sym('x0', nx);
% w = {w{:}, Xk};
% lbw = [lbw; x0_lb];
% ubw = [ubw; x0_ub];
% w0 = [w0; xg];
% 
% t = t0;
% dt = (tf-t0)/N;
% % Formulate the NLP
% for k=0:N-1
%     % New NLP variable for the control
%     Uk = casadi.MX.sym(['U_' num2str(k)], nu);
%     w = {w{:}, Uk};
%     lbw = [lbw; u_lb];
%     ubw = [ubw; u_ub];
%     w0 = [w0;  ug];
% 
%     % Integrate till the end of the interval
%     [xf, qf] = F(t, t0, tf, Xk, Uk, p);
%     Xk_end = xf;
%     J = J + qf;
%     
%     % next step
%     t = t + dt;
% 
%     % New NLP variable for state at end of interval
%     Xk = casadi.MX.sym(['X_' num2str(k+1)], nx);
%     w = {w{:}, Xk};
%     if k==N-1
%         lbw = [lbw; xf_lb];
%         ubw = [ubw; xf_ub];
%     else
%         lbw = [lbw; x_lb];
%         ubw = [ubw; x_ub];
%     end
%     w0 = [w0; xg];
% 
%     % Add equality constraint
%     g = [g, {Xk_end-Xk}];
%     lbg = [lbg; zeros(nx, 1)];
%     ubg = [ubg; zeros(nx, 1)];
% end
% 
% J = -Xk(2);
% 
% % Create an NLP solver
% prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
% solver = casadi.nlpsol('solver', 'ipopt', prob);
% 
% % Solve the NLP
% sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
% w_opt = full(sol.x);
% 
% % Plot the solution
% x1_opt = w_opt(3:4:end);
% x2_opt = w_opt(4:4:end);
% x3_opt = w_opt(5:4:end);
% u_opt = w_opt(6:4:end);
% tgrid = linspace(0, full(w_opt(2)), N+1);
% 
% figure(1)
% subplot(411); hold on;
% plot(tgrid, x1_opt)
% subplot(412); hold on;
% plot(tgrid, x2_opt)
% subplot(413); hold on;
% plot(tgrid, x2_opt)
% subplot(414); hold on;
% stairs(tgrid, [u_opt; nan])

