classdef ast_sign < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_sign(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = sign(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('sign(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end