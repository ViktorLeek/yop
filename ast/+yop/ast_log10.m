classdef ast_log10 < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_log10(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        function ast(obj)
            fprintf('log10(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end