classdef ast_log < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_log(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = log(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('log(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end