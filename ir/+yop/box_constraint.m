classdef box_constraint < yop.ast_node
    % just a ast_node to be able to draw. Should be changed!
    properties
        variable
        bound
    end
    methods
        function obj = box_constraint(variable, bound)
            obj.variable = variable;
            obj.bound = bound;
        end
        
        function draw(obj)
            fprintf([obj.name, '(variable, bound)\n']);
            
            begin_child(obj);
            draw(obj.variable);
            end_child(obj);
            
            last_child(obj);
            draw(obj.bound);
            end_child(obj);
        end
    end
end