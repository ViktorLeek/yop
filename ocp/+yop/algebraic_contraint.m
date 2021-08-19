classdef algebraic_contraint < yop.node
    properties
        expr
    end
    methods
        function obj = algebraic_contraint(expr)
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf('algebraic_contraint(expr)\n');
            
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
        
    end
end