classdef ast_ceil < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_ceil(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        function ast(obj)
            fprintf('ceil(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end