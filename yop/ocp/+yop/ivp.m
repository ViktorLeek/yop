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
            vars = yop.ocp.find_special_nodes(obj.equations{:});
            obj.classify_variables(vars);
            
            srf = yop.ocp.to_srf(obj.equations);
            [box, nbox] = yop.ocp.separate_box(srf);
            obj.set_box(box);
            
            [odes, algs, eqs, ieqs] = yop.ocp.sort_nonbox(nbox);
            obj.vectorize_ode(odes);
            obj.vectorize_alg(eqs);
            
            if ~isempty(algs) || ~isempty(ieqs)
                error(yop.msg.unexpected_error);
            end
            
            [odee, alge] = obj.set_ivp_functions();
            dae = struct;
            dae.t = mx_vec(obj.independent);
            dae.x = mx_vec(obj.states);
            dae.z = mx_vec(obj.algebraics);
            dae.p = mx_vec(obj.parameters);
            dae.ode = odee;
            dae.alg = dae.z+2.0;
            
            opts = struct;
            opts.output_t0 = true;
            opts.grid = linspace( ...
                obj.independent_initial.lb, ...
                obj.independent_final.lb, ...
                200);
            opts.print_stats = true;
            
            x0=[];
            for k=1:length(obj.states)
                x0 = [x0(:); obj.states(k).lb0(:)];
            end
            
            z0=[];
            for k=1:length(obj.algebraics)
                z0 = [z0(:); obj.algebraics(k).lb0(:)];
            end
            
            p=[];
            for k=1:length(obj.parameters)
                p = [p(:); obj.parameters(k).lb(:)];
            end
            
            F = casadi.integrator('F', 'idas', dae, opts);
            res = F('x0', x0, 'p', p);
            
        end
        
        function obj = vectorize_alg(obj, algs)
            alg_expr = [];
            for k=1:length(algs)
                alg_expr = [alg_expr; algs{k}.lhs(:)];
            end
            % g(t,x,z,p) == 0, should perharps be other way around, but it
            % does not matter, and it is consistent with canonlization of
            % OCPs.
            obj.alg = yop.ocp_rel(yop.ast_eq(alg_expr, 0));
        end
        
        function [ode_expr, alg_expr] = set_ivp_functions(obj)
            args = { ...
                mx_vec(obj.independent), ...
                mx_vec(obj.states), ...
                mx_vec(obj.algebraics), ...
                mx_vec(obj.parameters) ...
                };
            set_mx(obj.variables);
            ode_expr = fw_eval(obj.ode.ast.rhs);
            alg_expr = fw_eval(obj.alg.ast.lhs);
            obj.ode.fn = casadi.Function('ode', args, {ode_expr});
            obj.alg.fn = casadi.Function('ode', args, {alg_expr});
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
        
        function obj = set_box(obj, box)
            [box_t, box_t0, box_tf] = yop.ocp.timed_box(box);
            obj.parse_box(box_t, 'lb', 'ub');
            obj.parse_box(box_t0, 'lb0', 'ub0');
            obj.parse_box(box_tf, 'lbf', 'ubf');
            obj.set_box_bounds(); % Includes processing default values
        end
        
        function obj = parse_box(obj, box, lb, ub)
            for k=1:length(box)
                re = yop.reaching_elems(box{k}.lhs);
                bnd = yop.prop_num(box{k}.rhs);
                var = obj.find_variable(re.var.id);
                switch class(box{k})
                    case 'yop.ast_eq' % var == bnd
                        var.(ub)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                        var.(lb)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    case {'yop.ast_le', 'yop.ast_lt'} % var < bnd
                        var.(ub)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    case {'yop.ast_ge', 'yop.ast_gt'} % var > bnd
                        var.(lb)(re.reaching_idx) = ...
                            yop.get_subexpr(bnd, re.expr_elem);
                    otherwise
                        error('[Yop] Error: Wrong constraint class.');
                end
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
            vars = [obj.independent(:).', ...
                obj.independent_initial(:).', ...
                obj.independent_final(:).', ...
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
                obj.independent = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'independent variable']);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'initial bound for the independent variable']);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_var(t);
            else
                error(['[Yop] Error: An OCP can only have one ' ...
                    'terminal bound for the independent variable']);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ocp_var.empty(1,0);
            end
            obj.states(end+1) = yop.ocp_var(x);
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_var.empty(1,0);
            end
            obj.algebraics(end+1) = yop.ocp_var(z);
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_var.empty(1,0);
            end
            obj.parameters(end+1) = yop.ocp_var(p);
        end
        
        function obj = set_box_bounds(obj)
            obj.set_box_bnd('independent', 'lb', 'ub', ...
                yop.defaults().independent_lb, ...
                yop.defaults().independent_ub);
            
            obj.set_box_bnd('independent_initial', 'lb', 'ub', ...
                yop.defaults().independent_lb0, ...
                yop.defaults().independent_ub0);
            
            obj.set_box_bnd('independent_final', 'lb', 'ub', ...
                yop.defaults().independent_lbf, ...
                yop.defaults().independent_ubf);
            
            obj.set_box_bnd('states', 'lb', 'ub', ...
                yop.defaults().state_lb, yop.defaults().state_ub);
            
            obj.set_box_bnd('states', 'lb0', 'ub0', ...
                yop.defaults().state_lb0, yop.defaults().state_ub0);
            
            obj.set_box_bnd('states', 'lbf', 'ubf', ...
                yop.defaults().state_lbf, yop.defaults().state_ubf);
            
            obj.set_box_bnd('algebraics', 'lb', 'ub', ...
                yop.defaults().algebraic_lb, yop.defaults().algebraic_ub);
            
            obj.set_box_bnd('parameters', 'lb', 'ub', ...
                yop.defaults().parameter_lb, yop.defaults().parameter_ub);
        end 
        
        
        function obj = set_box_bnd(obj, var_field, lb_field, ...
                ub_field, lb_def, ub_def)
            for v=obj.(var_field)
                
                ub = v.(ub_field); % The full bound vector
                not_set = isnan(ub); % The elements not set
                
                bd = v.ub(not_set); % The ub for t = (t0, tf)
                bd(isnan(bd)) = ub_def; % Default for unsets
                
                ub(not_set) = bd; % Timed bound
                v.(ub_field) = ub; % Set variable timed bound
                
                % Same procedure for lb
                lb = v.(lb_field);
                not_set = isnan(lb);
                bd = v.lb(not_set);
                bd(isnan(bd)) = lb_def;
                lb(not_set) = bd;
                v.(lb_field) = lb;
            end
        end
        
    end
end