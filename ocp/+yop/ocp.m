classdef ocp < handle
    properties
        objective
        constraints
        variables
        box
        equality
        inequality
        independent
        independent_initial
        independent_final
        states
        algebraics
        controls
        parameters
        alg_eqs  % algebraic equations
        odes
    end
    methods
        function obj = ocp()
        end
        
        function obj = min(obj, objective)
            obj.objective = objective;
        end
        
        function obj = max(obj, objective)
            obj.objective = -objective;  % negate to get min problem
        end
        
        function obj = st(obj, varargin)
            obj.constraints = varargin;
        end
        
        function obj = build(obj)
            obj.parse_variables();
            obj.parse_constraints();
        end
        
        function obj = parse_variables(obj)
            vars = yop.get_variables({obj.objective, obj.constraints{:}});
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
                        
                    case 'yop.ast_algebraic'
                        obj.add_algebraic(vk);
                        
                    case 'yop.ast_control'
                        obj.add_control(vk);
                        
                    case 'yop.ast_parameter'
                        obj.add_parameter(vk);
                        
                    otherwise
                        error('[Yop] Error: Unknown variable type');
                        
                end
            end
        end
        
        function obj = parse_constraints(obj)            
            cons = yop.to_srf(obj.constraints);
            
            for k=1:length(cons)
                ck = yop.classify_constraint(cons{k});
                switch class(ck)
                    case 'yop.differential_constraint'
                        obj.add_ode(ck);
                        
                    case 'yop.algebraic_contraint'
                        obj.add_algebraic_eq(ck);
                        
                    case 'yop.box_initial_equality'
                        var = get_variable(ck.var);
                        
                        switch class(var)
                            case 'yop.ast_independent'
                                % What type of constraint is this??
                                % initial and ast_independent??
                                % Should be t(t0) - have to increase
                                % ocp_independent class.
                            case 'yop.ast_independent_initial'
                            case 'yop.ast_independent_final'
                            case 'yop.ast_state'
                            case 'yop.ast_algebraic'
                            case 'yop.ast_control'
                            case 'yop.ast_parameter'
                            otherwise
                                error('[Yop] Error: Unknown data type.');
                        end
                        
                    case 'yop.box_equality'
                        
                        % Need to match which variable the box constraint
                        % concernes. So the underlying variable of th
                        % contraint is found, and based on its type the
                        % same variable is found in ocp object. Once a
                        % match between the variables is found, it is
                        % possible to extract the box constraint values
                        % and set in the ocp_[type] object.
                        var = get_variable(ck.var);
                        
                        switch class(var)
                            case 'yop.ast_independent'
                            case 'yop.ast_independent_initial'
                                t = obj.get_matching_independent_initial(var);
                                t.set_ub(ck);
                                t.set_lb(ck);
                            case 'yop.ast_independent_final'
                            case 'yop.ast_state'
                            case 'yop.ast_algebraic'
                            case 'yop.ast_control'
                            case 'yop.ast_parameter'
                            otherwise
                                %error('[Yop] Error: Unknown data type.');
                        end
                        
                    case 'yop.box_final_equality'
                    case 'yop.box_initial_upper'
                    case 'yop.box_upper'
                    case 'yop.box_final_upper'
                    case 'yop.box_initial_lower'
                    case 'yop.box_lower'
                    case 'yop.box_final_lower'
                    case 'yop.inequality_constraint'
                    case 'yop.equality_constraint'
                    otherwise
                        error('[Yop] Error: Unknown constraint type.')
                end
            end
        end
        
        function v = get_matching_independent(obj, v2m)
            % v2m - variable to match
            for v = obj.independent
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_independent_initial(obj, v2m)
            % v2m - variable to match
            for v = obj.independent_initial
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_independent_final(obj, v2m)
            % v2m - variable to match
            for v = obj.independent_final
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_state(obj, v2m)
            % v2m - variable to match
            for v = obj.states
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_algebraic(obj, v2m)
            % v2m - variable to match
            for v = obj.algebraics
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_control(obj, v2m)
            % v2m - variable to match
            for v = obj.controls
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function v = get_matching_parameters(obj, v2m)
            % v2m - variable to match
            for v = obj.independent
                if isequal(v.var, v2m)
                    return;
                end
            end
            error('[Yop] Error: Could not match variable.')
        end
        
        function obj = add_independent(obj, t)
            if isempty(obj.independent)
                obj.independent = yop.ocp_independent(t);
            else
                obj.independent(end+1) = yop.ocp_independent(t);
            end
        end
        
        function obj = add_independent_initial(obj, t)
            if isempty(obj.independent_initial)
                obj.independent_initial = yop.ocp_independent_initial(t);
            else
                obj.independent_initial(end+1) = ...
                    yop.ocp_independent_initial(t);
            end
        end
        
        function obj = add_independent_final(obj, t)
            if isempty(obj.independent_final)
                obj.independent_final = yop.ocp_independent_final(t);
            else
                obj.independent_final(end+1) = ...
                    yop.ocp_independent_final(t);
            end
        end
        
        function obj = add_state(obj, x)
            if isempty(obj.states)
                obj.states = yop.ocp_state(x);
            else
                obj.states(end+1) = yop.ocp_state(x);
            end
        end
        
        function obj = add_algebraic(obj, z)
            if isempty(obj.algebraics)
                obj.algebraics = yop.ocp_algebraic(z);
            else
                obj.algebraics(end+1) = yop.ocp_algebraic(z);
            end
        end
        
        function obj = add_algebraic_eq(obj, eq)
            if isempty(obj.alg_eqs)
                obj.alg_eqs = yop.ocp_algebraic_eq(eq);
            else
                obj.alg_eqs(end+1) = yop.ocp_algebraic_eq(eq);
            end
        end
        
        function obj = add_control(obj, u)
            if isempty(obj.controls)
                obj.controls = yop.ocp_control(u);
            else
                obj.controls(end+1) = yop.ocp_control(u);
            end
        end
        
        function obj = add_parameter(obj, p)
            if isempty(obj.parameters)
                obj.parameters = yop.ocp_parameter(p);
            else
                obj.parameters(end+1) = yop.ocp_parameter(p);
            end
        end
        
        function obj = add_box(obj, bc)
            obj.box = [obj.box(:)', {bc}];
        end
        
        function obj = add_equality(obj, bc)
            obj.equality = [obj.equality(:)', {bc}];
        end
        
        function obj = add_inequality(obj, bc)
            obj.inequality = [obj.inequality(:)', {bc}];
        end
        
        function obj = add_ode(obj, ode)
            if isempty(obj.odes)
                obj.odes = yop.ocp_ode(ode.var, ode.expr);
            else
                obj.odes(end+1) = yop.ocp_ode(ode.var, ode.expr);
            end
        end
        
    end
end