classdef ast_sinh < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_sinh(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = sinh(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('sinh(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end