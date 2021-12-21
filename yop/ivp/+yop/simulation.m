classdef simulation < handle
    properties
        name
        
        % Variables
        independent 
        independent0
        independentf
        states
        algebraics
        parameters
        
        % dynamics
        ode
        alg
        
        % User inputs
        ode_eqs
        alg_eqs
    end
    
    properties (Hidden)
        visited % Visited nodes in the topological sort
    end
    
    methods
        function obj = simulation(varargin)
            obj.independent  = yop.ivp_var.empty(1,0);
            obj.independent0 = yop.ivp_var.empty(1,0);
            obj.independentf = yop.ivp_var.empty(1,0);
            obj.states       = yop.ivp_var.empty(1,0);
            obj.algebraics   = yop.ivp_var.empty(1,0);
            obj.parameters   = yop.ivp_var.empty(1,0);
            obj.add(varargin{:});
        end
        
        function obj = add(obj, varargin)
            for k=1:length(varargin)
                obj.parse_eq(varargin{k});
            end
        end
        
        function res = solve(obj, varargin)
            ip = inputParser();
            ip.FunctionName = "yop.simulation/solve";
            ip.addParameter('solver', yop.defaults.ivp_solver);
            ip.addParameter('opts', []);
            ip.parse(varargin{:});
            usr = ip.Results;
            
            obj.augment_system();
            obj.vectorize_dynamics();
            
            switch usr.solver
                case 'idas'
                    dae = struct( ...
                        't', obj.independent.mx_vec(), ...
                        'x', obj.states.mx_vec(), ...
                        'z', obj.algebraics.mx_vec(), ...
                        'p', obj.parameters.mx_vec(), ...
                        'ode', fw_eval(obj.ode.rhs), ...
                        'alg', fw_eval(obj.alg.rhs) ...
                        );
                    
                    if isempty(usr.opts)
                        opts = struct( ...
                            'output_t0', true, ...
                            'grid', linspace(obj.t0, obj.tf, yop.defaults.ivp_sol_points), ...
                            'print_stats', false ...
                            );
                    else
                        opts = usr.opts;
                    end
                    
                    F = casadi.integrator('F', 'idas', dae, opts);
                    res = F('x0', obj.x0, 'z0', obj.z0, 'p', obj.p0);
                    
                case 'cvodes'
                    if obj.is_dae()
                       error(yop.error.ode_solver_for_dae_problem());
                    end
                    
%                     dae = @(t,x) full([obj.)
                case 'collocation'
                case 'rk'
                case {'ode45', 'ode23', 'ode113', 'ode78', 'ode89'}
                case {'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'ode15i'}
                    
                    mx_fun = obj.ode_fn();
                    odefun = @(t,x) full(mx_fun(t,x,obj.p0));
                    options = odeset('Mass', obj.mass_matrix);
                    res = ode15s(odefun, ...
                        [obj.t0, obj.tf], ...
                        [obj.x0, obj.z0]', ...
                        options);
                    
                otherwise
                    error(yop.error.ivp_solver_not_recognized());
            end
            
        end
        
        function parse_eq(obj, eq)
            obj.find_variables(eq);
            ssr = yop.to_ssr(eq);
            for k=1:length(ssr)
                obj.classify_eq(ssr{k});
            end
        end
        
        function find_variables(obj, expression)
            % Find and add all special nodes of the expression
            [tsort, N, obj.visited] = ...
                topological_sort(expression, obj.visited);
            
            % Find special nodes, maintain topological order
            for n=1:N                
                switch class(tsort{n})
                    case 'yop.ast_independent'
                        obj.add_independent(tsort{n});
                    case 'yop.ast_independent_initial'
                        obj.add_independent0(tsort{n});
                    case 'yop.ast_independent_final'
                        obj.add_independentf(tsort{n});
                    case 'yop.ast_state'
                        obj.add_state(tsort{n});
                    case {'yop.ast_algebraic', 'yop.ast_control'}
                        obj.add_algebraic(tsort{n});
                    case 'yop.ast_parameter'
                        obj.add_parameter(tsort{n});
                    otherwise
                        continue
                end
            end
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ivp_var(t);
            else
                error(yop.error.multiple_independent_variables());
            end
        end
        
        function obj = add_independent0(obj, t)
            if isempty(obj.independent0)
                obj.independent0 = yop.ivp_var(t);
            else
                error(yop.error.multiple_independent_initial());
            end
        end
        
        function obj = add_independentf(obj, t)
            if isempty(obj.independentf)
                obj.independentf = yop.ivp_var(t);
            else
                error(yop.error.multiple_independent_final());
            end
        end
        
        function obj = add_state(obj, x)
            obj.states(end+1) = yop.ivp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            obj.algebraics(end+1) = yop.ivp_var(z);
        end
        
        function obj = add_parameter(obj, p)
            obj.parameters(end+1) = yop.ivp_var(p);
        end
        
        function classify_eq(obj, ssr)
            
            lhs = ssr.lhs;
            rhs = ssr.rhs;
            t0 = yop.initial_timepoint();
            
            isa_eq = isa(ssr, 'yop.ast_eq');
            
            lhs_num = isa_numeric(lhs);
            isa_der_lhs = isa_der(lhs);
            [~, id_lhs] = isa_variable(lhs);
            isa_state_lhs = isa_state(lhs);
            isa_control_lhs = isa_control(lhs);
            isa_algebraic_lhs = isa_algebraic(lhs);
            isa_parameter_lhs = isa_parameter(lhs);
            [lhs_istp, lhs_tp] = isa_timepoint(lhs);
            isa_independent0_lhs = isa_independent0(lhs);
            isa_independentf_lhs = isa_independentf(lhs);
            
            rhs_num = isa_numeric(rhs);
            isa_der_rhs = isa_der(rhs);
            [~, id_rhs] = isa_variable(rhs);
            isa_state_rhs = isa_state(rhs);
            isa_control_rhs = isa_control(rhs);
            isa_algebraic_rhs = isa_algebraic(rhs);
            isa_parameter_rhs = isa_parameter(rhs);
            [rhs_istp, rhs_tp] = isa_timepoint(rhs);
            isa_independent0_rhs = isa_independent0(rhs);
            isa_independentf_rhs = isa_independentf(rhs);
                
            if isa_der_lhs && isa_state_lhs && isa_eq
                % der(x) == expr
                obj.ode_eqs{end+1} = ssr;
                
            elseif isa_der_rhs && isa_state_rhs && isa_eq
                % expr == der(x)
                c = get_constructor(ssr);
                obj.ode_eqs{end+1} = c(ssr.rhs, ssr.lhs);
                
            elseif lhs_istp && lhs_tp==t0 && rhs_num && isa_eq && (isa_state_lhs || isa_control_lhs || isa_algebraic_lhs)
                % v(t0) == num
                var = obj.find_variable(id_lhs);
                var.iv = yop.prop_num(rhs);
                
            elseif lhs_num && rhs_istp && rhs_tp==t0 && isa_eq && (isa_state_rhs || isa_control_rhs || isa_algebraic_rhs)
                % num == v(t0)
                var = obj.find_variable(id_rhs);
                var.iv = yop.prop_num(lhs);
                
            elseif isa_independent0_lhs && rhs_num && isa_eq
                % t0 == num
                var = obj.find_variable(id_lhs);
                var.iv = yop.prop_num(rhs);
                
            elseif isa_independent0_rhs && lhs_num && isa_eq
                % num == t0
                var = obj.find_variable(id_rhs);
                var.iv = yop.prop_num(lhs);
                
            elseif isa_independentf_lhs && rhs_num && isa_eq
                % tf == num
                var = obj.find_variable(id_lhs);
                var.iv = yop.prop_num(rhs);
                
            elseif isa_independentf_rhs && lhs_num && isa_eq
                % num == tf
                var = obj.find_variable(id_rhs);
                var.iv = yop.prop_num(lhs);
                
            elseif isa_parameter_lhs && rhs_num && isa_eq
                % p == num
                var = obj.find_variable(id_lhs);
                var.iv = yop.prop_num(rhs);
                
            elseif lhs_num && isa_parameter_rhs && isa_eq
                % num == p
                var = obj.find_variable(id_rhs);
                var.iv = yop.prop_num(lhs);
                
            elseif isa_eq
                % Algebraic equation
                obj.alg_eqs{end+1} = canonicalize(ssr).lhs;
                
            else
                % error
                error(yop.error.failed_to_parse_ivp_eq());
                
            end
        end
        
        function ocp_var = find_variable(obj, id)
            for ocp_var = obj.variables()
                if ocp_var.ast.id == id
                    return;
                end
            end
            error(yop.error.failed_to_find_variable(id));
        end
        
        function vars = variables(obj)
            vars = [ ...
                obj.independent0(:).', ...
                obj.independentf(:).', ...
                obj.independent(:).', ...
                obj.states(:).', ...
                obj.algebraics(:).', ...
                obj.parameters(:).' ...
                ];
        end
        
        function obj = vectorize_dynamics(obj)
            obj.sort_states();
            
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
        
        function ids = sort_states(obj)
            % Change state order so that they come in id order
            [ids, idx] = sort(obj.get_state_ids());
            obj.states = obj.states(idx);
        end
        
        function ids = get_state_ids(obj)
            nx = length(obj.states);
            ids = zeros(nx,1);
            for k=1:nx
                ids(k) = obj.states(k).ast.id;
            end
        end
        
        function obj = set_dynamics_fn(obj)
            args = obj.mx_args();
            obj.variables.set_mx();
            obj.ode.fn = ...
                casadi.Function('ode', args, {fw_eval(obj.ode.rhs)});
            obj.alg.fn = ...
                casadi.Function('alg', args, {fw_eval(obj.alg.rhs)});
        end
        
        function args = mx_args(obj)
            args = { ...                
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.parameters) ...
                };
        end
        
        function args = ode_args(obj)
            args = { ...
                mx_vec(obj.independent), ...
                [mx_vec(obj.states); mx_vec(obj.algebraics)], ...
                mx_vec(obj.parameters) ...
                };
        end
        
        function rhs = ode_fn(obj)
            args = obj.ode_args();
            obj.variables.set_mx();
            expr = [fw_eval(obj.ode.rhs); fw_eval(obj.alg.rhs)];
            rhs = casadi.Function('ode', args, {expr});
        end
        
        function augment_system(obj)
            
            if isempty(obj.independent)
                obj.add_independent(yop.independent());
            end
            
            if isempty(obj.independent0)
                error(yop.error.start_of_integration_missing());
            end
            
            if isempty(obj.independentf)
                error(yop.error.end_of_integration_missing());
            end
            
        end
        
        function iv = t0(obj)
            if isempty(obj.independent0.iv)
                error(yop.error.initial_value_missing(obj.independent0));
            else
                iv = obj.independent0.iv;
            end
        end
        
        function iv = tf(obj)
            if isempty(obj.independentf.iv)
                error(yop.error.final_value_missing(obj.independentf));
            else
                iv = obj.independentf.iv;
            end
        end
        
        function iv = x0(obj)
            iv = [];
            for x=obj.states
                if isempty(x.iv)
                    error(yop.error.initial_value_missing(x));
                else
                    iv(end+1) = x.iv;
                end
            end
        end
        
        function iv = z0(obj)
            iv = [];
            for z=obj.algebraics
                if isempty(z.iv)
                    iv(end+1) = yop.defaults.algebraic_guess;
                else
                    iv(end+1) = z.iv;
                end
            end
        end
        
        function iv = p0(obj)
            % Constant value
            iv = [];
            for p=obj.parameters
                if isempty(p.iv)
                    error(yop.error.initial_value_missing(p));
                else
                    iv(end+1) = p.iv;
                end
            end
        end
        
        function n = n_x(obj)
            n = length(obj.states);
        end 
        
        function n = n_z(obj)
            n = length(obj.algebraics);
        end 
        
        function bool = is_ode(obj)
            bool = ~isempty(obj.ode) && isempty(obj.alg);
        end
        
        function bool = is_dae(obj)
            bool = ~isempty(obj.ode) && ~isempty(obj.alg);
        end
        
        function M = mass_matrix(obj)
            M = diag([ones(1,length(obj.states)), ...
                zeros(1,length(obj.algebraics))]);
        end
        
    end
end