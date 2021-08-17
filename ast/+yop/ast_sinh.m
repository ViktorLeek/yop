classdef ast_sinh < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_sinh(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = sinh(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('sinh(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end