classdef ast_erfinv < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_erfinv(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = erfin(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('erfinv(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end