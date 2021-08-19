classdef ocp < handle
    properties
        objective
        constraints
        variables = {}
        box = {}
        equality = {}
        inequality = {}
        independent = {}
        independent_initial = {}
        independent_final = {}
        states = {}
        algebraics = {}
        controls = {}
        parameters = {}
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
                    case 'yop.box_initial_equality'
                    case 'yop.box_equality'
                    case 'yop.box_final_equality'
                    case 'yop.box_initial_upper'
                    case 'yop.box_upper'
                    case 'yop.box_final_upper'
                    case 'yop.box_initial_lower'
                    case 'yop.box_lower'
                    case 'yop.box_final_lower'
                    case 'yop.inequality_constraint'
                    case 'yop.equality_contraint'
                    otherwise
                        error('[yop] Error: Unknown constraint type.')
                end
            end
        end        
        
        function obj = add_independent(obj, t)
            obj.independent = {obj.independent{:}, t};
        end
        
        function obj = add_independent_initial(obj, t)
            obj.independent_initial = {obj.independent_initial{:}, t};
        end
        
        function obj = add_independent_final(obj, t)
            obj.independent_final = {obj.independent_final{:}, t};
        end
        
        function obj = add_state(obj, x)
            obj.states = {obj.states{:}, x};
        end
        
        function obj = add_algebraic(obj, z)
            obj.algebraics = {obj.algebraics{:}, z};
        end
        
        function obj = add_control(obj, u)
            obj.controls = {obj.controls{:}, u};
        end
        
        function obj = add_parameter(obj, p)
            obj.parameters = {obj.parameters{:}, p};
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
            
            for k=1:length(obj.states)
                sk = obj.states{k};
                if sk.var == ode.var
                    sk.ode = ode.expr;
                    
                elseif 1==1
                    % der(x(s)) = expr
                    
                end
            end
        end
        
    end
end