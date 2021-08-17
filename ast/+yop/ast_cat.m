classdef ast_cat < yop.ast_expression
    properties
        d
        args
    end
    methods
        function obj = ast_cat(d, varargin)
            obj@yop.ast_expression();
            obj.d = d;
            obj.args = varargin;
            tmp = cell(size(varargin));
            for k=1:length(varargin)
                tmp{k} = ones(size(varargin{k}));
            end
            obj.dim = size(cat(d, tmp{:}));
        end
        
        function value = evaluate(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = evaluate(obj.args{k});
            end
            value = cat(evaluate(obj.d), tmp{:});
        end
        
        function v = forward(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = value(obj.args{k});
            end
            obj.m_value = cat(value(obj.d), tmp{:});
            v = obj.m_value;
        end
        
        function draw(obj)
            % every vararg is enumerated: "a1, a2, ..., aN, "
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['cat(dim, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            draw(obj.d);
            end_child(obj);
            for k=1:(length(obj.args)-1)
                begin_child(obj);
                draw(obj.args{k});
                end_child(obj);
            end
            last_child(obj);
            draw(obj.args{end});
            end_child(obj);
        end
        
        function [topsort, visited] = topological_sort(obj, topsort, visited)
            % Topological sort of expression graph by a dfs.
            
            % Initialize if second and third args are empty
            if nargin == 1
                topsort = {};
                visited = [];
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited] = topological_sort(obj.d, topsort, visited);
            
            for k=1:length(obj.args)
                [topsort, visited] = topological_sort(...
                    obj.args{k}, ...
                    topsort, ...
                    visited ...
                    );
            end
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
        
        
    end
end