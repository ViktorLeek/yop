classdef ast_mod < yop.ast_expression
    properties
        a
        m
    end
    methods
        function obj = ast_mod(a, m)
            obj@yop.ast_expression();
            obj.a = a;
            obj.m = m;
            obj.dim = size(mod(ones(size(a)), ones(size(m))));
        end
        
        function value = evaluate(obj)
            value = mod(evaluate(obj.a), evaluate(obj.m));
        end
        
        function v = forward(obj)
            obj.m_value = mod(value(obj.a), value(obj.m));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('mod(a, m)\n');
            
            begin_child(obj);
            draw(obj.a);
            end_child(obj);
            
            last_child(obj);
            draw(obj.m);
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
                topological_sort(obj.a, topsort, visited, n_elem);
            
            [topsort, visited, n_elem] = ...
                topological_sort(obj.m, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end