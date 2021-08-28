classdef ast_functioncall < yop.ast_node
    properties
        args
        dim
    end
    methods
        function obj = ast_functioncall(varargin)
            obj.args = varargin;
        end
    end
    methods % printing
        function ast(obj)
            % Print name with variable number of arguments
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf([obj.name, '(', str(1:end-2), ')\n']);
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