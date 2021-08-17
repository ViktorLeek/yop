classdef ast_norm < yop.ast_expression
    properties
        expr
        p
        nargs
    end
    methods
        function obj = ast_norm(expr, p)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.nargs = nargin;
            obj.dim = [1, 1];
            if nargin==2
                obj.p = p;
            end
        end
        
        function value = evaluate(obj)
            if obj.nargs == 1
                value = norm(evaluate(obj.expr));
            else
                value = norm(evaluate(obj.expr), evaluate(obj.p));
            end
        end
        
        function draw(obj)
            if obj.nargs == 1
                fprintf('norm(expr)\n');
                last_child(obj);
                draw(obj.expr);
                end_child(obj);
            else
                fprintf('norm(expr, p)\n');
                begin_child(obj);
                draw(obj.expr);
                end_child(obj);
                last_child(obj);
                draw(obj.p);
                end_child(obj);
            end
        end
    end
end