function nlp = dms_tmp(ocp, N, rk_steps)

if nargin == 1
    rk_steps = yop.defaults.rk4_steps;
end

%% OCP Variables
t = casadi.MX.sym('t');
t0 = casadi.MX.sym('t0');
tf = casadi.MX.sym('tf');

x_vec = [];
x_cell = cell(length(ocp.states), 1);
for k=1:length(ocp.states)
    xk = casadi.MX.sym( ...
        ocp.states(k).var.name, ...
        size(ocp.states(k).var, 1), ...
        size(ocp.states(k).var, 2));
    x_cell{k} = xk;
    x_vec = [x_vec(:); xk(:).'];
end
nx = numel(x_vec);

u_vec = [];
u_cell = cell(length(ocp.controls), 1);
for k=1:length(ocp.controls)
    uk = casadi.MX.sym( ...
        ocp.controls(k).var.name, ...
        size(ocp.controls(k).var, 1), ...
        size(ocp.controls(k).var, 2));
    u_cell{k} = uk;
    u_vec = [u_vec(:); uk(:).'];
end
nu = numel(u_vec);

p_vec = [];
p_cell = cell(length(ocp.parameters), 1);
for k=1:length(ocp.parameters)
    pk = casadi.MX.sym( ...
        ocp.parameters(k).var.name, ...
        size(ocp.parameters(k).var, 1), ...
        size(ocp.parameters(k).var, 2));
    p_cell{k} = pk;
    p_vec = [p_vec(:); pk(:).'];
end
np = numel(p_vec);



%% ODE
% Because the order in which the states come in the ocp variable vector and
% in the ode formulation might be different, it is necessary to derive a
% function for the ode lhs as well (for the variables that is). For
% instance, if the ocp state vector is x=[x1;x2] and the  the dynamics 
% are formualted as der([x2; x1]) == expr, then it necessary to permute the
% state vector, so that the constraint is parameterized correctly.
ode_var = []; ode_expr = [];
ocp.set_variables(t, t0, tf, x_cell, u_cell, p_cell);
for k=1:length(ocp.odes)
    lhs = evaluate(ocp.odes(k).var);
    rhs = forward_evaluate(ocp.odes(k).expr);
    ode_var = [ode_var(:); lhs(:).'];
    ode_expr = [ode_expr(:); rhs(:).'];
end
ocp.reset_variables();
ode_lhs = casadi.Function('ode_lhs', {x_vec}, {ode_var});
ode_rhs = casadi.Function('ode_rhs', {t, x_vec, u_vec, p_vec}, {ode_expr});

% Integrator for the dynamics
F = yop.rk4_integrator(ode_rhs, size(x_vec), size(u_vec), size(p_vec), rk_steps);


%% NLP variables
nlp = struct;

nlp.t = cell(N+1, 1);
nlp.t0 = t0;
nlp.tf = tf;
dt = (tf-t0)/N; % Grid step
nlp.t{1} = t0;
nlp.t{end} = tf;
for k=2:N
    nlp.t{k} = nlp.t{k-1} + dt;
end

nlp.x = cell(N+1,1);
for k=1:(N+1)
    nlp.x{k} = casadi.MX.sym(['x_' num2str(k)], nx);
end

nlp.u = cell(N, 1);
for k=1:N
    nlp.u{k} = casadi.MX.sym(['u_' num2str(k)], nu);
end

nlp.p = p_vec;

%% box constraints
[t0_lb, t0_ub, tf_lb, tf_ub] = independent_bound(ocp);
nlp.t0_lb = t0_lb;
nlp.t0_ub = t0_ub;
nlp.tf_lb = tf_lb;
nlp.tf_ub = tf_ub;

[x0_lb, x0_ub, x_lb, x_ub, xf_lb, xf_ub] = state_bound(ocp);
nlp.x_lb = repmat(x_lb, N+1, 1);
nlp.x_ub = repmat(x_ub, N+1, 1);
nlp.x_lb(1:nx) = x0_lb;
nlp.x_ub(1:nx) = x0_ub;
nlp.x_lb(end-nx+1:end) = xf_lb;
nlp.x_ub(end-nx+1:end) = xf_ub;

[u0_lb, u0_ub, u_lb, u_ub, uf_lb, uf_ub] = control_bound(ocp);
nlp.u_lb = repmat(u_lb, N, 1);
nlp.u_ub = repmat(u_ub, N, 1);
nlp.u_lb(1:nu) = u0_lb;
nlp.u_ub(1:nu) = u0_ub;
nlp.u_lb(end-nu+1:end) = uf_lb;
nlp.u_ub(end-nu+1:end) = uf_ub;

[p_lb, p_ub] = parameter_bound(ocp);
nlp.p_lb = p_lb;
nlp.p_ub = p_ub;

%% Objective function
J = parameterize_expression(ocp.objective);


%% Integrate trajectory
nlp.xf = cell(N, 1);
for k=1:N
    nlp.xf{k} = F(nlp.t{k}, nlp.t{k+1}, nlp.x{k}, nlp.u{k}, nlp.p);
end

%% Continuity/defect constraints
% ode_lhs should be a function in x only. However, 
nlp.def = cell(N,1);
for k=1:N
    nlp.def{k} = ode_lhs(nlp.x{k+1})-nlp.xf{k};
end

%% NLP Objective
nlp.J = 0; % This part should be integrated
Jfn = casadi.Function('J', {t, x_vec, u_vec, p_vec}, {evaluate(ocp.objective)});
nlp.J = Jfn(nlp.t{end}, nlp.x{end}, nlp.u{end}, nlp.p);

% Compile
nlp.w = vertcat(nlp.t0, nlp.tf, nlp.x{:}, nlp.u{:}, nlp.p);
nlp.w_lb = vertcat(nlp.t0_lb, nlp.tf_lb, nlp.x_lb, nlp.u_lb, nlp.p_lb);
nlp.w_ub = vertcat(nlp.t0_ub, nlp.tf_ub, nlp.x_ub, nlp.u_ub, nlp.p_ub);
nlp.g = vertcat(nlp.def{:});

nlp.ode_lhs = ode_lhs;
nlp.ode_rhs = ode_rhs;
nlp.Jfn = Jfn;
nlp.F = F;
% Restore values

end

function value = parameterize_expression(expr)

[topsort, n_elem] = topological_sort(obj.objective);

% Parameterize timepoints first, and after do not call forward on them as
% it can alter the value.
idx = [];
for k=1:n_elem
    if isa(topsort{k}, 'yop.ast_timepoint')
        % Set variables to their timepoint value
        % Evaluate expression
        % Set node value
        fn = obj.ode.fn(expr);
        [t,x,u,p] = 
        topsort{k}.m_value = 
    else
        idx(end+1) = k;
    end
end

% Integrate all integrations. This comes second as they might contain
% timepoints
idx2 = [];
for k=idx
    if isa(topsort{k}, 'yop.ast_int')
        % Integrera
    else
        idx2(end+1) = k;
    end
end

% Set the problem variables to the general OCP ones that can be used in a
% function generator.
ocp.set_variables(t,t0,tf,x,u,p);

% Finally, evaluate the expression in order to obatin an expression that
% can be used in a function generator with parameter list: t,x,z,u,p,w,
% where w is the full nlp vector.
for k=idx2
    forward(topsort{k});
end

value = topsort{k}.value;

ocp.reset_variables();

end






























% t = casadi.MX.sym('t');
% for k=1:length(ocp.independent)
%     ocp.independent(k).store_value();    
%     ocp.independent(k).set_value(t);
% end
% 
% t0 = casadi.MX.sym('t0');
% for k=1:length(ocp.independent_initial)
%     ocp.independent_initial(k).store_value();    
%     ocp.independent_initial(k).set_value(t0);
% end
% 
% tf = casadi.MX.sym('tf');
% for k=1:length(ocp.independent)
%     ocp.independent_final(k).store_value();    
%     ocp.independent_final(k).set_value(tf);
% end
% 
% x_vec = [];
% x_cell = cell(length(ocp.states), 1);
% for k=1:length(ocp.states)
%     ocp.states(k).store_value();
%     tmp = ocp.states(k).var;
%     xk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
%     x_cell{k} = xk;
%     ocp.states(k).set_value(xk);
%     x_vec = [x_vec(:); xk(:).'];
% end
% nx = numel(x_vec);
% 
% u_vec = [];
% u_cell = cell(length(ocp.controls), 1);
% for k=1:length(ocp.controls)
%     ocp.controls(k).store_value();
%     tmp = ocp.controls(k).var;
%     uk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
%     u_cell{k} = uk;
%     ocp.controls(k).set_value(uk);
%     u_vec = [u_vec(:); uk(:).'];
% end
% nu = numel(u_vec);
% 
% p_vec = [];
% p_cell = cell(length(ocp.parameters), 1);
% for k=1:length(ocp.parameters)
%     ocp.parameters(k).store_value();
%     tmp = ocp.parameters(k).var;
%     pk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
%     p_cell{k} = pk;
%     ocp.parameters(k).set_value(pk);
%     p_vec = [p_vec(:); pk(:).'];
% end
% np = numel(p_vec);

% exec = cell(n_elem, 1);
% cnt = 0;
% J_int = {};
% J_tp = {};
% for k=1:n_elem
%     switch class(topsort{k})
%         case 'yop.ast_int'
%             J_int = {J_int{:}, topsort{k}};
%         case 'yop.ast_timepoint'
%             J_tp = {J_tp{:}, topsort{k}};
%         otherwise
%             cnt = cnt + 1;
%             exec{cnt} = topsort{k};
%     end
% end
% exec = exec{1:cnt}; 
