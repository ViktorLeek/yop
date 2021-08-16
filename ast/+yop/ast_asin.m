classdef ast_asin < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_asin(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = asin(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('asin(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end