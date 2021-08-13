classdef ast_ctranspose < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_ctranspose(expr)
            obj.expr = expr;
            obj.dim = size(ctranspose(ones(size(expr))));
        end
        
        function value = evaluate(obj)
            value = ctranspose(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('ctranspose(obj)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
        
    end
end