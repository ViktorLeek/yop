classdef ast_acos < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_acos(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = acos(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('acos(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end