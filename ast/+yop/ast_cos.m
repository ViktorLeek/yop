classdef ast_cos < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_cos(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = cos(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('cos(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end