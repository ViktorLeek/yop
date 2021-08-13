classdef ast_floor < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_floor(expr)
            % Only accepts floor with single argument input
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            value = floor(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('floor(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end