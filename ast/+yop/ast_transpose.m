classdef ast_transpose < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_transpose(expr)
            obj.expr = expr;
            obj.dim = size(transpose(ones(size(expr))));
        end
        
        function value = evaluate(obj)
            value = transpose(evaluate(obj.expr));
        end
    
        function ast(obj)
            fprintf('transpose(obj)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end