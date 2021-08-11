classdef ast_reshape < yop.ast_node
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
    end
    methods % printing
        function ast(obj)
            str = [];
            for k=1:length(obj.szs)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['reshape(expr, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.expr);
            end_child(obj);
            for k=1:(length(obj.szs)-1)
                begin_child(obj);
                ast(obj.szs{k});
                end_child(obj);
            end
            last_child(obj);
            ast(obj.szs{end});
            end_child(obj);
        end
    end
end