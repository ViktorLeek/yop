classdef ast_reshape < yop.ast_expression
    properties
        expr
        szs
    end
    methods
        function obj = ast_reshape(expr, varargin)
            obj@yop.ast_expression();
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
        
        function v = forward(obj)
            tmp = cell(size(obj.szs));
            for k=1:length(tmp)
                tmp{k} = value(obj.szs{k});
            end
            obj.m_value = reshape(value(obj.expr), tmp{:});
            v = obj.m_value;
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
        
        function [topsort, visited, n_elem] = ...
                topological_sort(obj, topsort, visited, n_elem)
            % Topological sort of expression graph by a dfs.
            
            if nargin == 1
                % Start new sort
                visited = [];
                topsort = cell( ...
                    yop.constants().topsort_preallocation_size, 1);
                n_elem = 0;
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited, n_elem] = ...
                topological_sort(obj.expr, topsort, visited, n_elem);
            
            for k=1:length(obj.args)
                % probably unnecessary, as szs are expected to be numerics
                % but could change in the future.
                [topsort, visited, n_elem] = topological_sort( ...
                    obj.szs{k}, ...
                    topsort, ...
                    visited, ...
                    n_elem ...
                    );
            end
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end