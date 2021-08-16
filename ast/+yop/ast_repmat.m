classdef ast_repmat < yop.ast_expression
    properties
        expr
        args
    end
    methods
        function obj = ast_repmat(expr, varargin)
            obj.expr = expr;
            obj.args = varargin;
            obj.dim = repmat(ones(size(expr)), varargin{:});
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.args{k});
            end
            value = repmat(evaluate(obj.expr), tmp{:});
        end
        
        function ast(obj)
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['repmat(A, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.expr);
            end_child(obj);
            for k=1:(length(obj.args)-1)
                begin_child(obj);
                ast(obj.args{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.args{end});
            end_child(obj);
        end
    end
end