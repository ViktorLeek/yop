classdef ast_cosh < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_cosh(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = cosh(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('cosh(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end