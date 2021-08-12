classdef ast_atan < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_atan(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        function ast(obj)
            fprintf('atan(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end