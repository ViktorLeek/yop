classdef differential_constraint < yop.node
    properties
        var
        expr
    end
    methods
        function obj = differential_constraint(var, expr)
            obj.var = var;
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf('differential_constraint(var, expr)\n');
            
            begin_child(obj);
            draw(obj.var);
            end_child(obj);
            
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
        
    end
end