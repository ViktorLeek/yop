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
            obj.constraints = yop.to_srf(obj.constraints);
            
            % Classify constraints
            constraints = yop.to_srf(varargin);
            for k=1:length(constraints)
                [type, con] = get_constraint(constraints{k});
                switch type
                    case 'box'
                        obj.add_box(con);
                        
                    case 'equality'
                        obj.add_equality(con);
                        
                    case 'inequality'
                        obj.add_inequality(con);
                        
                    case 'differential'
                        obj.add_differential(con);
                        
                    case 'algebraic'
                        obj.add_algebraic(con);
                        
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