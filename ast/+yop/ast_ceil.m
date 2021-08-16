classdef ast_ceil < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_ceil(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = ceil(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('ceil(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end