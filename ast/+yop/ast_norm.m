classdef ast_norm < yop.ast_node
    properties
        expr
        p
        n_args
    end
    methods
        function obj = ast_norm(expr, p)
            obj.expr = expr;
            obj.n_args = 1;
            obj.dim = [1, 1];
            if nargin==2
                obj.p = p;
                obj.n_args = 2;
            end
        end
        function ast(obj)
            if obj.n_args == 1
                fprintf('norm(expr)\n');
                last_child(obj);
                ast(obj.expr);
                end_child(obj);
            else
                fprintf('norm(expr, p)\n');
                begin_child(obj);
                ast(obj.expr);
                end_child(obj);
                last_child(obj);
                ast(obj.p);
                end_child(obj);
            end
        end
    end
end