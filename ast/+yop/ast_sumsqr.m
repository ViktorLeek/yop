classdef ast_sumsqr < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_sumsqr(expr)
            obj.expr = expr;
            % dim is over simplified. To properly determine size it is
            % necessary to inspect the number of outputs the user expects.
            % Here however, the casadi approach of returning a single
            % variable is taken.
            % see: "doc casadi.GenericMatrixCommon.sumsqr"
            obj.dim = [1, 1];
        end
        
        function value = evaluate(obj)
            value = sumsqr(evaluate(obj.expr));
        end
        
        function ast(obj)
            fprintf('sumsqr(expr)\n');
            last_child(obj);
            ast(obj.expr);
            end_child(obj);
        end
    end
end