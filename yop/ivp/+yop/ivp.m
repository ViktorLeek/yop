classdef ivp < handle
    properties
        equations
        independent
        independent_initial
        independent_final
        states
        algebraics
        parameters
        
        ode
        alg
    end
    methods
        
        function obj = ivp(varargin)
            obj.equations = varargin;
        end
        
        function sol = solve(obj)
            vars = yop.ivp.find_variables(obj.equations{:});
            obj.classify_variables(vars);
            
            srf = yop.ocp.to_srf(obj.equations);
            [T, iv, dae_eq] = yop.ivp.split_iv(srf);
            obj.set_horizon(T);
            obj.set_iv(iv);
            
            [odes, algs] = yop.ivp.sort_eqs(dae_eq);
            obj.vectorize_ode(odes);
            obj.vectorize_alg(algs);
            
            [odee, alge] = obj.set_ivp_functions();
            dae = struct;
            dae.t = mx_vec(obj.independent);
            dae.x = mx_vec(obj.states);
            dae.z = mx_vec(obj.algebraics);
            dae.p = mx_vec(obj.parameters);
            dae.ode = odee;
            dae.alg = alge;
            
            d = 9;
            %             N = ceil(obj.tf-obj.t0)*200;
            N = ceil(obj.tf-obj.t0)*1;
            cp = 'legendre';
            [grid, tau] = obj.grid(N, d, cp);
            
            opts = struct;
            opts.output_t0 = true;
            opts.grid = grid;
            opts.print_stats = false;
            opts.reltol = 1e-3;
            
            x0=[];
            for k=1:length(obj.states)
                x0 = [x0(:); obj.states(k).lb(:)];
            end
            
            z0=[];
            for k=1:length(obj.algebraics)
                z0 = [z0(:); obj.algebraics(k).lb(:)];
            end
            
            p=[];
            for k=1:length(obj.parameters)
                p = [p(:); obj.parameters(k).lb(:)];
            end
            
            F = casadi.integrator('F', 'idas', dae, opts);
            res = F('x0', x0, 'p', p);
            
            sol = yop.ivp_sol( ...
                obj.variables, ...
                obj.mx_args, ...
                grid(1), ...
                grid(end), ...
                grid, ...
                full(res.xf), ...
                full(res.zf), ...
                p, N, d, cp ...
                );
        end
        
        function [g, tau] = grid(obj, N, d, cp)
            tau = full([0, casadi.collocation_points(d, cp)]);
            dt = (obj.tf-obj.t0)/N;
            g = [];
            tn = obj.t0;
            for n=1:N
                g = [g, tn + tau*dt];
                tn = tn + dt;
            end
            g(end+1) = obj.tf;
        end
        
        function obj = vectorize_alg(obj, algs)
            alg_expr = [];
            for k=1:length(algs)
                alg_expr = [alg_expr; algs{k}.lhs(:)];
            end
            % g(t,x,z,p) == 0
            obj.alg = yop.ocp_rel(yop.ast_eq(alg_expr, 0));
        end
        
        function [ode_expr, alg_expr] = set_ivp_functions(obj)
            set_mx(obj.variables);
            ode_expr = fw_eval(obj.ode.ast.rhs);
            alg_expr = fw_eval(obj.alg.ast.lhs);
            obj.ode.fn = casadi.Function('ode', obj.mx_args, {ode_expr});
            obj.alg.fn = casadi.Function('alg', obj.mx_args, {alg_expr});
        end
        
        function args = mx_args(obj)
            args = { ...
                mx_vec(obj.independent_initial), ...
                mx_vec(obj.independent_final), ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.parameters) ...
                };
        end
        
        function obj = classify_variables(obj, vars)
            for k=1:length(vars)
                vk = vars{k};
                switch class(vk)
                    case 'yop.ast_independent'
                        obj.add_independent(vk);
                    case 'yop.ast_independent_initial'
                        obj.add_independent_initial(vk);
                    case 'yop.ast_independent_final'
                        obj.add_independent_final(vk);
                    case 'yop.ast_state'
                        obj.add_state(vk);
                    case {'yop.ast_algebraic', 'yop.ast_control'}
                        obj.add_algebraic(vk);
                    case 'yop.ast_parameter'
                        obj.add_parameter(vk);
                    otherwise
                        error('[Yop] Error: Unknown variable type');
                end
            end
            if isempty(obj.independent)
                obj.add_independent(yop.ast_independent('t'));
            end
            if isempty(obj.independent_initial)
                obj.add_independent_initial( ...
                    yop.ast_independent_initial('t0'));
            end
            if isempty(obj.independent_final)
                obj.add_independent_final( ...
                    yop.ast_independent_final('tf'));
            end
            if isempty(obj.states)
                obj.add_state(yop.ast_state('x', 0, 0));
            end
            if isempty(obj.algebraics)
                obj.add_algebraic(yop.ast_algebraic('z', 0, 0));
            end
            if isempty(obj.parameters)
                obj.add_parameter(yop.ast_parameter('p', 0, 0));
            end
        end
        
        function [T0, Tf] = set_horizon(obj, T)
            % T0 <= t <= TF
            for k=1:length(T)
                re = yop.reaching_elems(T{k}.lhs);
                bnd = yop.prop_num(T{k}.rhs);
                var = obj.find_variable(re.var.id);
                switch class(T{k})
                    case 'yop.ast_eq' % var == bnd
                        value = yop.get_subexpr(bnd, re.expr_elem);
                        var.ub(re.reaching_idx) = value;
                        var.lb(re.reaching_idx) = value;
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.ub(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.lb(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
            end
            [T0, Tf] = horizon(obj);
        end
        
        function [T0, Tf] = horizon(obj)
            T0 = t0(obj);
            Tf = tf(obj);
        end
        
        function t = t0(obj)
            if isnan(obj.independent.lb)
                t = obj.independent_initial.lb;
                assert(t==obj.independent_initial.ub, ...
                    yop.msg.ivp_ub_differ);
            else
                t = obj.independent.lb;
                assert(~isnan(t), yop.msg.ivp_no_start_time);
            end
        end
        
        function t = tf(obj)
            if isnan(obj.independent.ub)
                t = obj.independent_final.lb;
                assert(t==obj.independent_final.ub, ...
                    yop.msg.ivp_ub_differ);
            else
                t = obj.independent.ub;
                assert(~isnan(t), yop.msg.ivp_no_start_time);
            end
        end
        
        function obj = set_iv(obj, iv)
            for k=1:length(iv)
                re = yop.reaching_elems(iv{k}.lhs);
                bnd = yop.prop_num(iv{k}.rhs);
                var = obj.find_variable(re.var.id);
                switch class(iv{k})
                    case 'yop.ast_eq'
                        value = yop.get_subexpr(bnd, re.expr_elem);
                        var.ub(re.reaching_idx) = value;
                        var.lb(re.reaching_idx) = value;
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
                
%                 [bool, tp] = iv{k}.lhs.isa_timepoint();
%                 if bool % Should be scalar! otherwise error is good!
%                    if tp ~= t0(obj) && tp ~= yop.initial_timepoint
%                        error(yop.msg.ivp_t0_err);
%                    end
%                 end
            end
        end
        
        function ocp_var = find_variable(obj, id)
            for ocp_var = obj.variables
                if ocp_var.var.id == id
                    return;
                end
            end
            error('[Yop] Error: ID not found');
        end
        
        function vars = variables(obj)
            vars = [ ...
                obj.independent_initial(:).', ...
                obj.independent_final(:).', ...
                obj.independent(:).', ...
                obj.states(:).', ...
                obj.algebraics(:).', ...
                obj.parameters(:).'];
        end
        
        function obj = vectorize_ode(obj, odes)
            ode_var = []; ode_expr = [];
            for k=1:length(odes)
                ode_var = [ode_var(:); odes{k}.lhs(:)];
                ode_expr = [ode_expr(:); odes{k}.rhs(:)];
            end

            [re, nr] = obj.reaching_states(ode_var);
            
            for k=1:length(re)
                if numel(re(k).enum) ~= numel(re(k).reaching)
                    [~, idx] = setdiff(re(k).enum, re(k).reaching);
                    vk = re(k).var(idx);
                    ode_var = [ode_var(:); vk(:)];
                    ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                    
                    if yop.settings.warnings
                        warning(yop.msg.default_der(vk.name));
                    end
                end
            end
            
            for k=1:length(nr)
                vk = nr(k).var;
                ode_var = [ode_var(:); vk(:)];
                ode_expr = [ode_expr(:); zeros(size(vk(:)))];
                
                if yop.settings.warnings
                    warning(yop.msg.default_der(vk.name));
                end
            end
            
            i2e = obj.compute_permutation_vector(ode_var);
            
            ast = yop.ast_eq(ode_var(i2e), ode_expr(i2e));
            obj.ode = yop.ocp_rel(ast);
        end
        
        function [i2e, e2i] = compute_permutation_vector(obj, ode_var)
            [~, ~, output_idx] = obj.reaching_states(ode_var);
            [~, idx] = sort(output_idx);
            e2i = output_idx;
            i2e = idx;
        end
        
        function [re, nr, reaching_enumeration] = reaching_states(obj,expr)
            % Cell array with all states (ocp_variable objects)
            xs = arrayfun(@(e) e.var, obj.states, 'UniformOutput', false);
            re(length(xs)) = yop.re_data();
            
            % Enumerate all variable elements
            e0 = 1; % Enumeration start value.
            for k=1:length(re)
                re(k).var = xs{k};
                [e0, idx] = re(k).enumerate(e0);
                obj.states(k).enum = idx;
            end
            
            % Propagate elements through the expression
            reaching_enumeration = propagate_value(expr);
            
            % Filter
            reaching_enumeration(~(isa_der(expr) & isa_state(expr))) = 0;
            
            % Store results
            for k=1:length(re)
                re(k).set_expr_elements_reached(reaching_enumeration);
            end
            
            % Split reaching and not reaching
            reaching = arrayfun(@(v) ~isempty(v.reaching), re);
            nr = re(~reaching);
            re = re( reaching);
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ivp_var(t);
            else
                error(['[Yop] Error: An IVP can only have one ' ...
                    'independent variable']);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ivp_var(t);
            else
                error(['[Yop] Error: An IVP can only have one ' ...
                    'initial bound for the independent variable']);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ivp_var(t);
            else
                error(['[Yop] Error: An IVP can only have one ' ...
                    'terminal bound for the independent variable']);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ivp_var.empty(1,0);
            end
            obj.states(end+1) = yop.ivp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ivp_var.empty(1,0);
            end
            obj.algebraics(end+1) = yop.ivp_var(z);
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ivp_var.empty(1,0);
            end
            obj.parameters(end+1) = yop.ivp_var(p);
        end
        
    end
    
    methods (Static)
        
        
        function vars = find_variables(varargin)
            vars = {};
            visited = [];
            for k=1:length(varargin)
                [tsort, N, visited]=topological_sort(varargin{k}, visited);
                for n=1:N
                    if isa(tsort{n}, 'yop.ast_variable')
                        vars{end+1} = tsort{n};
                    end
                end
            end
        end
        
        function [T, iv, rem] = split_iv(srf)
            % Split into time horizon equations, initial values and 
            % problem eqs.
            T = {};
            iv = {};
            rem = {};
            for n=1:length(srf)
                var_num = yop.ivp.isa_T(srf{n}.lhs, srf{n}.rhs);                
                num_var = yop.ivp.isa_T(srf{n}.rhs, srf{n}.lhs);
                isT = var_num | num_var;
                T{end+1} = yop.get_subrel(srf{n}, isT);
                
                var_num = yop.ivp.isa_iv(srf{n}.lhs, srf{n}.rhs);                
                num_var = yop.ivp.isa_iv(srf{n}.rhs, srf{n}.lhs);
                isiv = var_num | num_var;
                iv{end+1} = yop.get_subrel(srf{n}, isiv);
                
                rem{end+1} = yop.get_subrel(srf{n}, ~isiv & ~isT);
            end
            
            T = T(~cellfun('isempty', T));
            for k=1:length(T)
                T{k} = canonicalize_box(T{k});
            end
            T = yop.ocp.unique_box(T);
            
            iv = iv(~cellfun('isempty', iv));
            for k=1:length(iv)
                iv{k} = canonicalize_box(iv{k});
            end
            iv = yop.ocp.unique_box(iv);
            
            rem = rem(~cellfun('isempty', rem));
        end
        
        function boolv = isa_T(var_cand, num_cand)
            boolv = isa_independent(var_cand)  & isa_numeric(num_cand);
        end
        
        function boolv = isa_iv(var_cand, num_cand)
            boolv = (...
                (isa_state(var_cand)     & isa_timepoint(var_cand)) | ...
                (isa_algebraic(var_cand) & isa_timepoint(var_cand)) | ...
                isa_parameter(var_cand) ...
                ) & ~isa_der(var_cand) & ...
                isa_numeric(num_cand);
        end
        
        function [ode, alg] = sort_eqs(nbox)
            ode={}; 
            alg={};
            for k=1:length(nbox)
                if isa(nbox{k}, 'yop.ast_eq')
                    [~, ode_k, alg_k] = isa_ode(nbox{k});
                    ode{end+1} = ode_k;
                    alg{end+1} = alg_k; % is canonicalized
                    
                else 
                    error(yop.msg.ivp_relation);
                    
                end
            end
            ode = ode(~cellfun('isempty', ode));
            alg = alg(~cellfun('isempty', alg));
        end 
    end
end




































