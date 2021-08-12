classdef ast_sum < yop.ast_node
    properties
        expr
        args
    end
    methods
        function obj = ast_sum(expr, varargin)
            % ast_sum
            
            % This only works for cases where the optional arguments are
            % not ast_node
            obj.expr = expr;
            obj.args = varargin;  % May be changed depending on case
            obj.dim = size(sum(ones(size(expr)), varargin{:}));
        end
        function ast(obj)
            
            if isempty(obj.args)
                fprintf('sum(expr)\n');
                last_child(obj);
                ast(obj.expr);
                end_child(obj);
            else
                str = [];
                for k=1:length(obj.args)
                    str = [str, 'a', num2str(k), ', '];
                end
                
                fprintf(['sum(expr, ', str(1:end-2), ')\n']);
                
                % print expr
                begin_child(obj);
                ast(obj.expr);
                end_child(obj);
                
                % print the rest of the arguments
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
end