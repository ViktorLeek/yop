classdef ast_erf < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_erf(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = erf(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('erf(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end