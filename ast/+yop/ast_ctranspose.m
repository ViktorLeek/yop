classdef ast_ctranspose < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_ctranspose(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(ctranspose(ones(size(expr))));
        end
        
        function value = evaluate(obj)
            value = ctranspose(evaluate(obj.expr));
        end
        
        function draw(obj)
            fprintf('ctranspose(obj)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
        
    end
end