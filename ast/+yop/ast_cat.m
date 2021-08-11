classdef ast_cat < yop.ast_node
    properties
        d
        args
    end
    methods
        function obj = ast_cat(d, varargin)
            obj.d = d;
            obj.args = varargin;
            tmp = cell(size(varargin));
            for k=1:length(varargin)
                tmp{k} = ones(size(varargin{k}));
            end
            obj.dim = size(cat(d, tmp{:}));
        end
    end
    methods % printing
        function ast(obj)
            % every vararg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['cat(dim, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            ast(obj.d);
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