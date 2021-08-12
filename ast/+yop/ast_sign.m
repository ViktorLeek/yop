classdef ast_sign < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_sign(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        function ast(obj)
            fprintf('sign(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end