classdef ast_trace < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_trace(expr)
            obj.expr = expr;
            obj.dim = size(trace(ones(size(expr))));
        end
        function ast(obj)
            fprintf('trace(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end