classdef ast_trace < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_trace(expr)
            obj.expr = expr;
            obj.dim = size(trace(ones(size(expr))));
        end
        
        function value = evaluate(obj)
            value = trace(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('trace(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end