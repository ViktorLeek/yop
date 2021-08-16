classdef ast_reshape < yop.ast_expression
    properties
        expr
        szs
    end
    methods
        function obj = ast_reshape(expr, varargin)
            obj.expr = expr;
            obj.szs = varargin;
            obj.dim = size(reshape(ones(size(expr)), varargin{:}));
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.szs));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.szs{k});
            end
            value = reshape(evaluate(obj.expr), tmp{:});
        end
        
        function draw(obj)
            str = [];
            for k=1:length(obj.szs)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['reshape(expr, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            draw(obj.expr);
            end_child(obj);
            for k=1:(length(obj.szs)-1)
                begin_child(obj);
                draw(obj.szs{k});
                end_child(obj);
            end
            last_child(obj);
            draw(obj.szs{end});
            end_child(obj);
        end
    end
end