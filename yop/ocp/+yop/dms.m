function nlp = dms(ocp, N, rk_steps)

if nargin == 1
    rk_steps = yop.defaults.rk4_steps;
end

t = casadi.MX.sym('t');
for k=1:length(ocp.independent)
    ocp.independent(k).store_value();    
    ocp.independent(k).set_value(t);
end

t0 = casadi.MX.sym('t0');
for k=1:length(ocp.independent_initial)
    ocp.independent_initial(k).store_value();    
    ocp.independent_initial(k).set_value(t0);
end

tf = casadi.MX.sym('tf');
for k=1:length(ocp.independent)
    ocp.independent_final(k).store_value();    
    ocp.independent_final(k).set_value(tf);
end

x = [];
for k=1:length(ocp.states)
    ocp.states(k).store_value();
    tmp = ocp.states(k).var;
    xk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
    ocp.states(k).set_value(xk);
    x = [x(:); xk(:).'];
end
nx = numel(x);

u = [];
for k=1:length(ocp.controls)
    ocp.controls(k).store_value();
    tmp = ocp.controls(k).var;
    uk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
    ocp.controls(k).set_value(uk);
    u = [u(:); uk(:).'];
end
nu = numel(u);

p = [];
for k=1:length(ocp.parameters)
    ocp.parameters(k).store_value();
    tmp = ocp.parameters(k).var;
    pk = casadi.MX.sym(tmp.name, size(tmp, 1), size(tmp, 2));
    ocp.parameters(k).set_value(pk);
    p = [p(:); pk(:).'];
end
np = numel(p);

ode_var = []; ode_expr = [];
for k=1:length(ocp.odes)
    lhs = evaluate(ocp.odes(k).var);
    rhs = forward_evaluate(ocp.odes(k).expr);
    ode_var = [ode_var(:); lhs(:).'];
    ode_expr = [ode_expr(:); rhs(:).'];
end

Jfn = casadi.Function('J', {t, x, u, p}, {evaluate(ocp.objective)});
ode_lhs = casadi.Function('ode_lhs', {x}, {ode_var});
ode_rhs = casadi.Function('ode_rhs', {t, x, u, p}, {ode_expr});
F = yop.rk4_integrator(ode_rhs, size(x), size(u), size(p), rk_steps);

% NLP variables
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

nlp.p = p;

% box constraints
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

% Integrate
nlp.xf = cell(N, 1);
for k=1:N
    nlp.xf{k} = F(nlp.t{k}, nlp.t{k+1}, nlp.x{k}, nlp.u{k}, nlp.p);
end

% Continuity
% ode_lhs should be a function in x only. However, 
nlp.def = cell(N,1);
for k=1:N
    nlp.def{k} = ode_lhs(nlp.x{k+1})-nlp.xf{k};
end

% Objective
nlp.J = 0; % This part should be integrated
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


