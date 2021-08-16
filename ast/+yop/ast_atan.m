classdef ast_atan < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_atan(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = atan(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('atan(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end