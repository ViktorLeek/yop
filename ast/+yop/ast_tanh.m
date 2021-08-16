classdef ast_tanh < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_tanh(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = tanh(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('tanh(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end