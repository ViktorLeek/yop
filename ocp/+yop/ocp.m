classdef ocp < handle
    properties
        objective
        constraints
        variables = {}
        box = {}
        equality = {}
        inequality = {}
        differential = {}
        algebraic = {}
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
        
        function obj = get_variables(obj)
            obj.variables = ...
                yop.get_variables({obj.objective, obj.constraints{:}});
        end
        
        function obj = classify_constraints(obj)            
            cons = yop.to_srf(obj.constraints);
            
            for k=1:length(cons)
                ck = yop.classify_constraint(cons{k});
                switch class(ck)
                    case 'yop.differential_contraint'
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
        
        function obj = add_box(obj, bc)
            obj.box = [obj.box(:)', {bc}];
        end
        
        function obj = add_equality(obj, bc)
            obj.equality = [obj.equality(:)', {bc}];
        end
        
        function obj = add_inequality(obj, bc)
            obj.inequality = [obj.inequality(:)', {bc}];
        end
        
        function obj = add_differential(obj, bc)
            obj.differential = [obj.differential(:)', {bc}];
        end
        
        function obj = add_algebraic(obj, bc)
            obj.algebraic = [obj.algebraic(:)', {bc}];
        end
    end
end