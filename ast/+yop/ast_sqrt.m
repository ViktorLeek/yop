classdef ast_sqrt < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_sqrt(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = sqrt(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('sqrt(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end