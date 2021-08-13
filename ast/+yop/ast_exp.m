classdef ast_exp < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_exp(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = exp(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('exp(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end