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
        
        function v = forward(obj)
            tmp = cell(size(obj.args));
            for k=1:length(tmp)
                tmp{k} = value(obj.args{k});
            end
            obj.m_value = repmat(value(obj.expr), tmp{:});
            v = obj.m_value;
        end
        
        function draw(obj)
            str = [];
            for k=1:length(obj.args)
                str = [str, 'a', num2str(k), ', '];
            end
            fprintf(['repmat(A, ', str(1:end-2), ')\n']);
            
            begin_child(obj);
            draw(obj.expr);
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
            [topsort, visited]=topological_sort(obj.expr, topsort, visited);
            for k=1:length(obj.args)
                % probably unnecessary, as args are expected to be numerics
                % but could change in the future.
                [topsort, visited] = topological_sort( ...
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