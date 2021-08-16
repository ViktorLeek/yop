classdef ast_tan < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_tan(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = tan(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('tan(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end