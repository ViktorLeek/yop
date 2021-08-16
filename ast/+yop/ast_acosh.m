classdef ast_acosh < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_acosh(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = acosh(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('acosh(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end