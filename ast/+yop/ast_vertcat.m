classdef ast_vertcat < yop.ast_node
    properties
        args
    end
    methods
        function obj = ast_vertcat(varargin)
            obj.args = varargin;
            tmp = cell(size(varargin));
            for k=1:length(varargin)
                tmp{k} = ones(size(varargin{k}));
            end
            obj.dim = size(vertcat(tmp{:}));
        end
    end
    methods % printing
        function ast(obj)
            % every arg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['vertcat(', str(1:end-2), ')\n']);
            
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