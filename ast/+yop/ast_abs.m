classdef ast_abs < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_abs(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = abs(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('abs(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end