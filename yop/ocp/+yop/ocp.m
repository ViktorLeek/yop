classdef ocp < handle
    properties
        name
        
        % Variables
        independent
        independent0
        independentf
        states
        algebraics
        controls
        parameters
        
        % Objective
        objective
        
        % Dynamics
        ode
        alg
        
        % Constraints
        path
        path_hard
        path_ival
        point
        
        % User input
        ode_eqs
        alg_eqs
        ec_eqs
        ec_hard_eqs
        ec_ival_eqs
        ec_point_eqs
        iec_eqs
        iec_hard_eqs
        iec_ival_eqs
        iec_point_eqs
        
    end
    
    properties (Hidden) % internal properties
        snodes % Topological sort of special nodes
        tps  % Timepoints
        ints % Integrals
        ders % Derivativtes (time varying value)
        visited = [] % special nodes that has been visited - speedup
    end
    
    methods
        function obj = ocp(name)
            if nargin == 1
                obj.name = name;
            else
                obj.name = 'Optimal Control Problem';
            end
            obj.ode_eqs = {};
            obj.alg_eqs = {};
            obj.ec_eqs  = {};
            obj.ec_hard_eqs  = {};
            obj.ec_ival_eqs  = {};
            obj.ec_point_eqs = {};
            obj.iec_eqs = {};
            obj.iec_hard_eqs  = {};
            obj.iec_ival_eqs  = {};
            obj.iec_point_eqs = {};
            
            obj.snodes = yop.ocp_expr.empty(1,0);
            obj.tps  = yop.ocp_expr.empty(1,0);
            obj.ints = yop.ocp_expr.empty(1,0);
            obj.ders = yop.ocp_expr.empty(1,0);
            
            obj.independent  = yop.ocp_independent.empty(1,0);
            obj.independent0 = yop.ocp_independent0.empty(1,0);
            obj.independentf = yop.ocp_independentf.empty(1,0);
            obj.states       = yop.ocp_state_control.empty(1,0);
            obj.algebraics   = yop.ocp_algebraic.empty(1,0);
            obj.controls     = yop.ocp_state_control.empty(1,0);
            obj.parameters   = yop.ocp_parameter.empty(1,0);
            
        end
        
        function obj = min(obj, expr)
            obj.set_objective(expr);
        end
        
        function obj = max(obj, expr)
            obj.set_objective(-expr);
        end
        
        function obj = st(obj, varargin)
            for k=1:length(varargin)
                obj.parse_constraint(varargin{k});
            end
        end
        
        function obj = hard(obj, varargin)
            for k=1:length(varargin)
                varargin{k} = hard(varargin{k});
            end
            obj.st(varargin{:});
        end
        
        function sol = solve(obj, varargin)
            ip = inputParser();
            ip.FunctionName = "yop.ocp/solve";
            ip.addParameter('intervals', yop.defaults.control_invervals);
            ip.addParameter('degree', yop.defaults.polynomial_degree);
            ip.addParameter('points', yop.defaults.collocation_points);
            ip.addParameter('guess', []);
            ip.parse(varargin{:});
            opts = ip.Results;
            
            N  = opts.intervals;
            d  = opts.degree;
            cp = opts.points;
            guess = opts.guess;
            
            % The system is augmented before box bounds are set in order to
            % use the correct default value for the augemented variables.
            obj.augment_system();
            obj.sort_states();
            obj.set_box_bounds();
            obj.set_objective_fn();
            obj.vectorize_dynamics();
            obj.set_point_con();
            obj.set_path_con();
            obj.set_hard_path_con();
            obj.set_ival_path_con();
            obj.set_special_functions();
            
            nlp = yop.direct_collocation(obj, N, d, cp);
            
            nlp_opts = struct;
            nlp_opts.ipopt.acceptable_tol = 1e-6;
            solver = casadi.nlpsol('solver', 'ipopt', ...
                struct('f', nlp.J, 'x', nlp.w, 'g', nlp.g), nlp_opts);
            nlp_sol = solver( ...
                'x0', ones(size(nlp.w)), ...
                'lbx', nlp.w_lb, ...
                'ubx', nlp.w_ub, ...
                'ubg', nlp.g_ub, ...
                'lbg', nlp.g_lb ...
                );
            
            w_opt = struct;
            w_opt.t0 = nlp.t0(nlp_sol.x);
            w_opt.tf = nlp.tf(nlp_sol.x);
            w_opt.t = nlp.t(nlp_sol.x);
            w_opt.x = nlp.x(nlp_sol.x);
            w_opt.z = nlp.z(nlp_sol.x);
            w_opt.u = nlp.u(nlp_sol.x);
            w_opt.p = nlp.p(nlp_sol.x);
            
            sol = yop.ocp_sol( ...
                obj.independent0, ...
                obj.independentf, ...
                obj.independent, ...
                obj.states, ...
                obj.algebraics, ...
                obj.controls, ...
                obj.parameters, ...
                obj.mx_vars(), ...
                w_opt, N, d, cp);
        end
        
        function augment_system(obj)
            
            % Augment system based on control parametrization
            % Step 1: Account for all control inputs
            for uk=obj.controls
                du = uk.ast.der;
                while ~isempty(du)
                    obj.add_unique_control(du);
                    du = du.der;
                end
            end
            
            % Step 2: Promote integrated controls to states and add
            %         augmenting equations
            keep = [];
            for k=1:length(obj.controls)
                uk = obj.controls(k);
                if ~isempty(uk.ast.der)
                    obj.states(end+1) = uk;
                    obj.ode_eqs{end+1} = ode(der(uk.ast)==uk.ast.der);
                else
                    keep(end+1) = k;
                end
            end
            obj.controls = obj.controls(keep);

            % Fill in blanks
            if isempty(obj.independent)
                obj.add_independent(yop.independent());
            end
            
            if isempty(obj.independent0)
                obj.add_independent0(yop.independent0());
            end
            
            if isempty(obj.independentf)
                obj.add_independentf(yop.independentf());
            end
            
        end
        
        function set_box_bounds(obj)
            
            % Time0
            if isempty(obj.independent0.ub) && isempty(obj.independent0.lb)
                obj.independent0.ub = yop.defaults.independent0_ub;
                obj.independent0.lb = yop.defaults.independent0_lb;
            end
                
            if isempty(obj.independent0.ub)
                obj.independent0.ub = inf;
            end
                
            if isempty(obj.independent0.lb)
                obj.independent0.lb = -inf;
            end
            
            % Timef
            if isempty(obj.independentf.ub) && isempty(obj.independentf.lb)
                obj.independentf.ub = yop.defaults.independentf_ub;
                obj.independentf.lb = yop.defaults.independentf_lb; 
            end
            
            if isempty(obj.independentf.ub)
                obj.independentf.ub = inf;
            end
            
            if isempty(obj.independentf.lb)
                % User should introduce t0 <= tf
                warning(['[Yop] Final time is unbounded from below. ', ...
                    'If this is intentional, consider introducing ''t0 <= tf''. ', ...
                    'If this was unintentional, set a lower bound for tf ' ...
                    '(''tf >= value'') and consider introduction the above mentioned ' ...
                    'constraint if t0 and tf can overlap.']);
                obj.independentf.lb = -inf; 
            end
            
            % Time
            % Enables constraints such as 't > 0, t < 10'
            if ~isempty(obj.independent.lb)
                t_min = obj.independent.lb;
                obj.independent0.lb = max(obj.independent0.lb, t_min);
                obj.independent0.ub = max(obj.independent0.ub, t_min);
                obj.independentf.fb = max(obj.independentf.lb, t_min);
                obj.independentf.fb = max(obj.independentf.ub, t_min);
            end
            
            if ~isempty(obj.independent.lb)
                t_max = obj.independent.ub;
                obj.independent0.lb = min(obj.independent0.lb, t_max);
                obj.independent0.ub = min(obj.independent0.ub, t_max);
                obj.independentf.fb = min(obj.independentf.lb, t_max);
                obj.independentf.fb = min(obj.independentf.ub, t_max);
            end
            
            % State
            for x=obj.states
                if isempty(x.ub)
                    x.ub = yop.defaults.state_ub;
                end
                if isempty(x.lb)
                    x.lb = yop.defaults.state_lb;
                end
                if isempty(x.ub0)
                    x.ub0 = x.ub;
                end
                if isempty(x.lb0)
                    x.lb0 = x.lb; 
                end
                if isempty(x.ubf)
                    x.ubf = x.ub;
                end
                if isempty(x.lbf)
                    x.lbf = x.lb;
                end
            end
            
            % Algebraics 
            for z=obj.algebraics
                if isempty(z.ub)
                    z.ub = yop.defaults.algebraic_ub;
                end
                if isempty(z.lb)
                    z.lb = yop.defaults.algebraic_lb;
                end
            end
            
            % Controls
            for u=obj.controls
                if isempty(u.ub)
                    u.ub = yop.defaults.control_ub;
                end
                if isempty(u.lb)
                    u.lb = yop.defaults.control_lb;
                end
                if isempty(u.ub0)
                    u.ub0 = u.ub;
                end
                if isempty(u.lb0)
                    u.lb0 = u.lb; 
                end
                if isempty(u.ubf)
                    u.ubf = u.ub;
                end
                if isempty(u.lbf)
                    u.lbf = u.lb;
                end
            end
            
            % Parameters
            for p=obj.parameters
                if isempty(p.ub)
                    p.ub = yop.defaults.parameter_ub;
                end
                if isempty(p.lb)
                    p.lb = yop.defaults.parameter_lb;
                end
            end
        end
        
        function set_objective(obj, expr)
            % Topological sort of expression in order to find variables,
            % timepoints, integrals and derivatives.
            tp_int = obj.find_special_nodes(expr);
            
            % Error handling - Test if timedependent variables reach the
            % objective without being in a timepoint or integral. By
            % making a second topological sort, where timepoints and
            % integrals have been removed, it is tested whether time
            % dependent variabels are reached in the ast. If a time
            % dependent variable is inside a timepoint or integral
            % it will not be visited by the topological sort.
            if yop.settings.errors
                visited_ = yop.get_ids(tp_int);
                [sort, N] = topological_sort(expr, visited_);
                err_nodes = {};
                for n=1:N
                    time_varying = ...
                        isa(sort{n}, 'yop.ast_state') || ...
                        isa(sort{n}, 'yop.ast_control') || ...
                        isa(sort{n}, 'yop.ast_algebraic') || ...
                        isa(sort{n}, 'yop.ast_independent');
                    if time_varying
                        err_nodes{end+1} = sort{n};
                    end
                end
                
                if ~isempty(err_nodes)
                    error(yop.error.timevarying_objective(err_nodes));
                end
            end
            
            % Passed error check - assign value
            obj.objective.ast = expr;
        end
        
        function obj = set_objective_fn(obj)
            args = { ...
                mx_vec(obj.independent0), ...
                mx_vec(obj.independentf), ...
                mx_vec(obj.parameters), ...
                mx_vec(obj.tps), ...
                mx_vec(obj.ints), ...
                };
            obj.set_mx();
            obj.snodes.set_mx();
            obj.objective.fn = ...
                casadi.Function('fn', args, {fw_eval(obj.objective.ast)});
        end
        
        function obj = set_special_functions(obj)
            args = obj.mx_args();
            obj.set_mx();
            obj.snodes.set_mx();
            for sn = obj.snodes
                mx_expr = fw_eval(sn.ast.expr);
                sn.fn = casadi.Function('fn', args, {mx_expr(:)});
            end
        end
        
        function vectorize_dynamics(obj)
            % Vectorize the equation
            n_ode = length(obj.ode_eqs);
            tmp_lhs = cell(n_ode, 1);
            tmp_rhs = cell(n_ode, 1);
            for k=1:length(obj.ode_eqs)
                tmp_lhs{k} = obj.ode_eqs{k}.lhs;
                tmp_rhs{k} = obj.ode_eqs{k}.rhs;
            end
            ode_lhs = vertcat(tmp_lhs{:});
            %ode_rhs = vertcat(tmp_rhs{:});
            
            % Test if all states are bound to an ode
            [~, ode_ids] = isa_state(ode_lhs);
            [ode_ids, idx] = sort(ode_ids);
            x_ids = obj.get_state_ids();
            if ~isequal(x_ids, ode_ids)
                state_ast = {};
                for id = setdiff(ode_ids, x_ids)
                    state_ast{end+1} = obj.find_variable(id);
                end
                error(yop.error.missing_state_derivative(state_ast));
            end
            
            % Change order of equations so that state vector and ode
            % equation order match
            obj.ode.lhs = ode_lhs(idx);
            obj.ode.rhs = vertcat(tmp_rhs{idx}); %ode_rhs(idx);
            
            % Algebraic equation
            nz = length(obj.alg_eqs);
            alg_rhs = cell(nz,1);
            for k=1:nz
                z_k = obj.alg_eqs{k};
                alg_rhs{k} = z_k.rhs - z_k.lhs;
            end
            obj.alg.lhs = zeros(nz, 1);
            obj.alg.rhs = vertcat(alg_rhs{:});
            
            % Compute symbolic functions of the dynamics
            obj.set_dynamics_fn();
        end
        
        function obj = set_dynamics_fn(obj)
            args = obj.mx_args();
            obj.set_mx();
            obj.snodes.set_mx();
            obj.ode.fn = ...
                casadi.Function('ode', args, {fw_eval(obj.ode.rhs)});
            obj.alg.fn = ...
                casadi.Function('alg', args, {fw_eval(obj.alg.rhs)});
        end
        
        function set_path_con(obj)            
            obj.set_mx();
            obj.snodes.set_mx();
            pc_vec = vertcat(obj.ec_eqs{:}, obj.iec_eqs{:});
            fn = casadi.Function('eq', obj.mx_args(), {fw_eval(pc_vec)});
            
            n_eq  = length(obj.ec_eqs);
            n_ieq = length(obj.iec_eqs);
            obj.path.fn = fn;
            obj.path.ub = zeros(n_eq+n_ieq, 1);
            obj.path.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function set_hard_path_con(obj)
            obj.set_mx();
            obj.snodes.set_mx();
            pc_vec = vertcat(obj.ec_hard_eqs{:}, obj.iec_hard_eqs{:});
            fn = casadi.Function('eq', obj.mx_args(), {fw_eval(pc_vec)});
            
            n_eq  = length(obj.ec_hard_eqs);
            n_ieq = length(obj.iec_hard_eqs);
            obj.path_hard.fn = fn;
            obj.path_hard.ub = zeros(n_eq+n_ieq, 1);
            obj.path_hard.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function set_point_con(obj)
            obj.set_mx();
            obj.snodes.set_mx();
            pc_vec = vertcat(obj.ec_point_eqs{:}, obj.iec_point_eqs{:});
            fn = casadi.Function('eq', obj.mx_args(), {fw_eval(pc_vec)});
            
            n_eq = length(obj.ec_point_eqs);
            n_ieq = length(obj.iec_point_eqs);
            obj.point.fn = fn;
            obj.point.ub = zeros(n_eq+n_ieq, 1);
            obj.point.lb = [zeros(n_eq,1); -inf(n_ieq,1)];
        end
        
        function set_ival_path_con(obj)
            args = obj.mx_args();
            obj.set_mx();
            obj.snodes.set_mx();
            
            info = yop.ocp_ival.empty(1,0);
            for k=1:length(obj.ec_ival_eqs)
                ik = obj.ec_ival_eqs{k};
                [t0, tf] = get_ival(ik);
                info(k).t0 = t0;
                info(k).tf = tf;
                info(k).ast = ik;
                info(k).ub = 0;
                info(k).lb = 0;
                expr = fw_eval(ik);
                info(k).fn = casadi.Function('eq', args, {expr(:)});
            end
            
            for k=1:length(obj.iec_ival_eqs)
                ik = obj.iec_ival_eqs{k};
                [t0, tf] = get_ival(ik);
                info(k).t0 = t0;
                info(k).tf = tf;
                info(k).ast = ik;
                info(k).ub = 0;
                info(k).lb = -inf;
                expr = fw_eval(ik);
                info(k).fn = casadi.Function('eq', args, {expr(:)});
            end
            
            obj.path_ival = info;
        end
        
        function args = mx_vars(obj)
            args = { ...
                mx_vec(obj.independent0), ...
                mx_vec(obj.independentf), ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters) ...
                };
        end
        
        function args = mx_args(obj)
            args = { ...
                mx_vec(obj.independent0), ...
                mx_vec(obj.independentf), ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.controls), ...
                mx_vec(obj.parameters), ...
                mx_vec(obj.tps), ...
                mx_vec(obj.ints), ...
                mx_vec(obj.ders) ...
                };
        end
        
        function set_mx(obj)
            obj.independent0.set_mx();
            obj.independentf.set_mx();
            obj.independent.set_mx();
            obj.states.set_mx();
            obj.algebraics.set_mx();
            obj.controls.set_mx();
            obj.parameters.set_mx();
            obj.snodes.set_mx();
        end
        
        function ids = get_state_ids(obj)
            nx = length(obj.states);
            ids = zeros(nx,1);
            for k=1:nx
                ids(k) = obj.states(k).ast.id;
            end
        end
        
        function ids = sort_states(obj)
            % Change state order so that they come in id order
            [ids, idx] = sort(obj.get_state_ids());
            obj.states = obj.states(idx);
        end
        
        function parse_constraint(obj, c)
            obj.find_special_nodes(c);
            ssr = yop.to_ssr(c);
            for k=1:length(ssr)
                obj.classify_constraint(ssr{k});
            end
        end
        
        function classify_constraint(obj, ssr)
            
            lhs = ssr.lhs;
            rhs = ssr.rhs;
            t0 = yop.initial_timepoint();
            tf = yop.final_timepoint();
            
            isa_eq = isa(ssr, 'yop.ast_eq');
            isa_le = isa(ssr, 'yop.ast_le') || isa(ssr, 'yop.ast_lt');
            isa_ge = isa(ssr, 'yop.ast_ge') || isa(ssr, 'yop.ast_gt');
            invariant = is_transcription_invariant(ssr);
            hard = is_hard(ssr);
            
            lhs_num = isa_numeric(lhs);
            is_ival_lhs = is_ival(lhs);
            [isa_der_lhs, der_id_lhs] = isa_der(lhs);
            [isa_var_lhs, id_lhs] = isa_variable(lhs);
            isa_fnh_lhs = isa(ssr.lhs, 'function_handle');
            isa_state_lhs = isa_state(lhs);
            isa_control_lhs = isa_control(lhs);
            [lhs_istp, lhs_tp] = isa_timepoint(lhs);
            
            rhs_num = isa_numeric(rhs);
            is_ival_rhs = is_ival(rhs);
            [isa_der_rhs, der_id_rhs] = isa_der(rhs);
            [isa_var_rhs, id_rhs] = isa_variable(rhs);
            isa_fnh_rhs = isa(ssr.rhs, 'function_handle');
            isa_state_rhs = isa_state(rhs);
            isa_control_rhs = isa_control(rhs);
            [rhs_istp, rhs_tp] = isa_timepoint(rhs);
            
            if (is_ival_lhs || is_ival_rhs) && isa_eq
                % interval equality constraint
                obj.ec_ival_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif (is_ival(lhs) || is_ival(rhs))
                % interval inequality constraint
                obj.iec_ival_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif isa_der_lhs && isa_state_lhs && isa_eq
                % der(x) == expr
                obj.ode_eqs{end+1} = ssr;
                obj.remove_state_der(der_id_lhs);
                
            elseif isa_der_rhs && isa_state_rhs && isa_eq
                % expr == der(x)
                c = get_constructor(ssr);
                obj.ode_eqs{end+1} = c(ssr.rhs, ssr.lhs);
                obj.remove_state_der(der_id_rhs);
                
            elseif is_alg(ssr) && isa_eq
                % alg(expr1 == expr2)
                obj.alg_eqs{end+1} = ssr;
                
            elseif lhs_istp && lhs_tp==t0 && isa_var_lhs && rhs_num && isa_eq && (isa_state_lhs || isa_control_lhs)
                % v(t0) == num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ub0 = bnd;
                var.lb0 = bnd;
                
            elseif lhs_istp && lhs_tp==t0 && isa_var_lhs && rhs_num && isa_le && (isa_state_lhs || isa_control_lhs)
                % v(t0) <= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ub0 = bnd;
                
            elseif lhs_istp && lhs_tp==t0 && isa_var_lhs && rhs_num && isa_ge && (isa_state_lhs || isa_control_lhs)
                % v(t0) >= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==t0 && isa_var_rhs && isa_eq && (isa_state_rhs || isa_control_rhs)
                % num == v(t0)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub0 = bnd;
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==t0 && isa_var_rhs && isa_le && (isa_state_rhs || isa_control_rhs)
                % num <= v(t0)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==t0 && isa_var_rhs && isa_ge && (isa_state_rhs || isa_control_rhs)
                % num >= v(t0)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub0 = bnd;
                
            elseif lhs_istp && lhs_tp==tf && isa_var_lhs && rhs_num && isa_eq && (isa_state_lhs || isa_control_lhs)
                % v(tf) == num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ubf = bnd;
                var.lbf = bnd;
                
            elseif lhs_istp && lhs_tp==tf && isa_var_lhs && rhs_num && isa_le && (isa_state_lhs || isa_control_lhs)
                % v(tf) <= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ub0 = bnd;
                
            elseif lhs_istp && lhs_tp==tf && isa_var_lhs && rhs_num && isa_ge && (isa_state_lhs || isa_control_lhs)
                % v(tf) >= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==tf && isa_var_rhs && isa_eq && (isa_state_rhs || isa_control_rhs)
                % num == v(tf)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub0 = bnd;
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==tf && isa_var_rhs && isa_le && (isa_state_rhs || isa_control_rhs)
                % num <= v(tf)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.lb0 = bnd;
                
            elseif lhs_num && rhs_istp && rhs_tp==tf && isa_var_rhs && isa_ge && (isa_state_rhs || isa_control_rhs)
                % num >= v(tf)
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub0 = bnd;
                
            elseif isa_var_lhs && rhs_num && isa_eq
                % v == num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ub = bnd;
                var.lb = bnd;
                
            elseif isa_var_lhs && rhs_num && isa_le
                % v <= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.ub = bnd;
                
            elseif isa_var_lhs && rhs_num && isa_ge
                % v >= num
                var = obj.find_variable(id_lhs);
                bnd = yop.prop_num(rhs);
                var.lb = bnd;
                
            elseif lhs_num && isa_var_rhs && isa_eq
                % num == v
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub = bnd;
                var.lb = bnd;
                
            elseif lhs_num && isa_var_rhs && isa_le
                % num <= v
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.lb = bnd;
                
            elseif lhs_num && isa_var_rhs && isa_ge
                % num >= v
                var = obj.find_variable(id_rhs);
                bnd = yop.prop_num(lhs);
                var.ub = bnd;
                
            elseif isa_var_lhs && isa_fnh_rhs && isa_eq
                % v == @
                var = obj.find_variable(id_lhs);
                var.ub = rhs;
                var.lb = rhs;
                
            elseif isa_var_lhs && isa_fnh_rhs && isa_le
                % v <= @
                var = obj.find_variable(id_lhs);
                var.ub = rhs;
                
            elseif isa_var_lhs && isa_fnh_rhs && isa_ge
                % v >= @
                var = obj.find_variable(id_lhs);
                var.lb = rhs;
                
            elseif isa_fnh_lhs && isa_var_rhs && isa_eq
                % @ == v
                var = obj.find_variable(id_rhs);
                var.ub = lhs;
                var.lb = lhs;
                
            elseif isa_fnh_lhs && isa_var_rhs && isa_le
                % @ <= v
                var = obj.find_variable(id_rhs);
                var.lb = lhs;
                
            elseif isa_fnh_lhs && isa_var_rhs && isa_ge
                % @ >= v
                var = obj.find_variable(id_rhs);
                var.ub = lhs;
                
            elseif isa_eq && invariant
                % Transcription invariant equality constraint
                obj.ec_point_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif isa_eq && hard
                % Transcription invariant equality constraint
                obj.ec_hard_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif isa_eq
                % Equality constraint
                obj.ec_eqs{end+1} = canonicalize(ssr).lhs;
            
            elseif (isa_le || isa_ge) && invariant
                % Transcription invariant inequality constraint
                obj.iec_point_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif (isa_le || isa_ge) && hard
                % Hard inequality constraint
                obj.iec_hard_eqs{end+1} = canonicalize(ssr).lhs;
                
            elseif isa_le || isa_ge
                % Inequality constraint
                obj.iec_eqs{end+1} = canonicalize(ssr).lhs;
                
            else
                % error
                error(yop.error.unknown_constraint());
                
            end
        end
        
        function ocp_var = find_variable(obj, id)
            
            % Here Matlab's short circuit of logical expressions is used to
            % avoid evaluating the id comparison if the variable is empty
            if ~isempty(obj.independent) && obj.independent.ast.id == id
                ocp_var = obj.independent;
                return;
            end
            
            if ~isempty(obj.independent0) && obj.independent0.ast.id == id
                ocp_var = obj.independent0;
                return;
            end
            
            if ~isempty(obj.independentf) && obj.independentf.ast.id == id
                ocp_var = obj.independentf;
                return;
            end
            
            for ocp_var = [obj.states, obj.controls]
                if ocp_var.ast.id == id
                    return;
                end
            end
            
            for ocp_var = obj.algebraics
                if ocp_var.ast.id == id
                    return;
                end
            end
            
            for ocp_var = obj.parameters
                if ocp_var.ast.id == id
                    return;
                end
            end
            
            error(yop.error.failed_to_find_variable(id));
        end
        
        function remove_state_der(obj, id)
            % Remove the state derivative node if all elements covered by
            % the derivative are states.
            for k=1:length(obj.ders)
                [bool, id_k] = isa_der(obj.ders(k).ast);
                if all(bool) && all(id_k == id)
                    if all(isa_state(obj.ders(k).ast))
                        to_remove = obj.ders(k);
                        obj.ders = [obj.ders(1:k-1), obj.ders(k+1:end)]; 
                        % Also need to remove it from special nodes vector
                        for n=1:length(obj.snodes)
                            if obj.snodes(n) == to_remove
                                obj.snodes = [obj.snodes(1:n-1), ...
                                    obj.snodes(n+1:end)];
                                to_remove.set_mx(); % evaluation for debug
                                return;
                            end
                        end
                    end
                end
            end
        end
        
        function tp_int = find_special_nodes(obj, expression)
            % Find and add all special nodes of the expression
            
            % Used to find errors in the objective function
            tp_int = {};
            
            % Sort all nodes
            [tsort, N, obj.visited] = ...
                topological_sort(expression, obj.visited);
            
            % Find special nodes, maintain topological order
            for n=1:N
                tn = tsort{n};
                if isa(tn, 'yop.ast_variable')
                    obj.add_variable(tn);
                    
                elseif isa(tn, 'yop.ast_timepoint')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.tp);
                    obj.snodes(end+1) = sn;
                    obj.tps(end+1) = sn;
                    tp_int{end+1} = tn;
                    
                elseif isa(tn, 'yop.ast_int')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.int);
                    obj.snodes(end+1) = sn;
                    obj.ints(end+1) = sn;
                    tp_int{end+1} = tn;
                    
                elseif isa(tn, 'yop.ast_der')
                    sn = yop.ocp_expr(tn, yop.ocp_expr.der);
                    obj.snodes(end+1) = sn;
                    obj.ders(end+1) = sn;
                end
            end
        end
        
        function add_variable(obj, v)
            switch class(v)
                case 'yop.ast_independent'
                    obj.add_independent(v);
                case 'yop.ast_independent_initial'
                    obj.add_independent0(v);
                case 'yop.ast_independent_final'
                    obj.add_independentf(v);
                case 'yop.ast_state'
                    obj.add_state(v);
                case 'yop.ast_algebraic'
                    obj.add_algebraic(v);
                case 'yop.ast_control'
                    obj.add_control(v);
                case 'yop.ast_parameter'
                    obj.add_parameter(v);
            end
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ocp_independent(t);
            else
                yop.error.multiple_independent_variables();
            end
        end
        
        function obj = add_independent0(obj, t)
            if isempty(obj.independent0)
                obj.independent0 = yop.ocp_independent0(t);
            else
                error(yop.error.multiple_independent_initial());
            end
        end
        
        function obj = add_independentf(obj, t)
            if isempty(obj.independentf)
                obj.independentf = yop.ocp_independentf(t);
            else
                error(yop.error.multiple_independent_final());
            end
        end
        
        function obj = add_state(obj, x)
            obj.states(end+1) = yop.ocp_state_control(x);
        end
        
        function obj = add_algebraic(obj, z)
            obj.algebraics(end+1) = yop.ocp_algebraic(z);
        end
        
        function obj = add_control(obj, u)
            obj.controls(end+1) = yop.ocp_state_control(u);
        end
        
        function obj = add_unique_control(obj, u)
            for uk = obj.controls
                if uk.ast.id == u.id
                    return;
                end
            end
            obj.add_control(u);
        end
        
        function obj = add_parameter(obj, p)
            obj.parameters(end+1) = yop.ocp_parameter(p);
        end
        
        function n = n_x(obj)
            n = length(obj.states);
        end
        
        function n = n_z(obj)
            n = length(obj.algebraics);
        end
        
        function n = n_u(obj)
            n = length(obj.controls);
        end
        
        function n = n_p(obj)
            n = length(obj.parameters);
        end
        
        function n = n_tp(obj)
            n = n_elem(obj.tps);
        end
        
        function n = n_int(obj)
            n = n_elem(obj.ints);
        end
        
        function n = n_der(obj)
            n = n_elem(obj.ders);
        end
        
        function bd = t0_ub(obj, t)
            if isa(obj.independent0.ub, 'function_handle')
                bd = obj.independent0.ub(t);
            else
                bd = obj.independent0.ub;
            end
        end
        
        function bd = t0_lb(obj, t)
            if isa(obj.independent0.lb, 'function_handle')
                bd = obj.independent0.lb(t);
            else
                bd = obj.independent0.lb;
            end
        end
        
        function bd = tf_ub(obj, t)
            if isa(obj.independentf.ub, 'function_handle')
                bd = obj.independentf.ub(t);
            else
                bd = obj.independentf.ub;
            end
        end
        
        function bd = tf_lb(obj, t)
            if isa(obj.independentf.lb, 'function_handle')
                bd = obj.independentf.lb(t);
            else
                bd = obj.independentf.lb;
            end
        end
        
        function bd = x0_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ub0, 'function_handle')
                    bd(end+1) = v.ub0(t);
                else
                    bd(end+1) = v.ub0;
                end
            end
            bd = bd(:);
        end
        
        function bd = x0_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lb0, 'function_handle')
                    bd(end+1) = v.lb0(t);
                else
                    bd(end+1) = v.lb0;
                end
            end
            bd = bd(:);
        end
        
        function bd = x_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ub, 'function_handle')
                    bd(end+1) = v.ub(t);
                else
                    bd(end+1) = v.ub;
                end
            end
            bd = bd(:);
        end
        
        function bd = x_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lb, 'function_handle')
                    bd(end+1) = v.lb(t);
                else
                    bd(end+1) = v.lb;
                end
            end
            bd = bd(:);
        end
        
        function bd = xf_ub(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.ubf, 'function_handle')
                    bd(end+1) = v.ubf(t);
                else
                    bd(end+1) = v.ubf;
                end
            end
            bd = bd(:);
        end
        
        function bd = xf_lb(obj, t)
            bd = [];
            for v = obj.states
                if isa(v.lbf, 'function_handle')
                    bd(end+1) = v.lbf(t);
                else
                    bd(end+1) = v.lbf;
                end
            end
            bd = bd(:);
        end
        
        function bd = u0_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ub0, 'function_handle')
                    bd(end+1) = v.ub0(t);
                else
                    bd(end+1) = v.ub0;
                end
            end
            bd = bd(:);
        end
        
        function bd = u0_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lb0, 'function_handle')
                    bd(end+1) = v.lb0(t);
                else
                    bd(end+1) = v.lb0;
                end
            end
            bd = bd(:);
        end
        
        function bd = u_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ub, 'function_handle')
                    bd(end+1) = v.ub(t);
                else
                    bd(end+1) = v.ub;
                end
            end
            bd = bd(:);
        end
        
        function bd = u_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lb, 'function_handle')
                    bd(end+1) = v.lb(t);
                else
                    bd(end+1) = v.lb;
                end
            end
            bd = bd(:);
        end
        
        function bd = uf_ub(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.ubf, 'function_handle')
                    bd(end+1) = v.ubf(t);
                else
                    bd(end+1) = v.ubf;
                end
            end
            bd = bd(:);
        end
        
        function bd = uf_lb(obj, t)
            bd = [];
            for v = obj.controls
                if isa(v.lbf, 'function_handle')
                    bd(end+1) = v.lbf(t);
                else
                    bd(end+1) = v.lbf;
                end
            end
            bd = bd(:);
        end
        
        function bd = p_ub(obj, t)
            bd = [];
            if ~isempty(obj.parameters)
                if isa(obj.parameters.ub, 'function_handle')
                    bd = obj.parameters.ub(t);
                else
                    bd = obj.parameters.ub;
                end
            end
            bd = bd(:);
        end
        
        function bd = p_lb(obj, t)
            bd = [];
            if ~isempty(obj.parameters)
                if isa(obj.parameters.lb, 'function_handle')
                    bd = obj.parameters.lb(t);
                else
                    bd = obj.parameters.lb;
                end
            end
            bd = bd(:);
        end
        
        function bd = z_ub(obj, t)
            bd = [];
            if ~isempty(obj.algebraics)
                if isa(obj.algebraics.ub, 'function_handle')
                    bd = obj.algebraics.ub(t);
                else
                    bd = obj.algebraics.ub;
                end
            end
            bd = bd(:);
        end
        
        function bd = z_lb(obj, t)
            bd = [];
            if ~isempty(obj.algebraics)
                if isa(obj.algebraics.lb, 'function_handle')
                    bd = obj.algebraics.lb(t);
                else
                    bd = obj.algebraics.lb;
                end
            end
            bd = bd(:);
        end
        
        function [bool, t0, tf] = fixed_horizon(obj)
            t0_ub = obj.independent0.ub;
            t0_lb = obj.independent0.lb;
            tf_ub = obj.independentf.ub;
            tf_lb = obj.independentf.lb;
            
            bool = t0_lb == t0_ub && tf_lb == tf_ub && ...
                all(~isinf([t0_lb, t0_ub, tf_lb, tf_ub]));
            
            t0 = t0_lb;
            tf = tf_lb;
        end
        
        function bool = has_path(obj)
            bool = ~isempty(obj.path.ub);
        end
        
        function bool = has_hard_path(obj)
            bool = ~isempty(obj.path_hard.ub);
        end
        
    end
end