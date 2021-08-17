classdef ast_exp < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_exp(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = exp(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('exp(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end