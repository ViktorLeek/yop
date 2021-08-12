classdef ast_tanh < yop.ast_node
    properties
        expr
    end
    methods
        function obj = ast_tanh(expr)
            obj.expr = expr;
            obj.dim = size(expr);
        end
        function ast(obj)
            fprintf('tanh(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end