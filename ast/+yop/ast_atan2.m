classdef ast_atan2 < yop.ast_expression
    properties
        y
        x
    end
    methods
        function obj = ast_atan2(y, x)
            obj@yop.ast_expression();
            obj.y = y;
            obj.x = x;
            obj.dim = size(atan2(ones(size(y)), ones(size(x))));
        end
        
        function value = evaluate(obj)
            value = atan2(evaluate(obj.y), evaluate(obj.x));
        end
        
        function v = forward(obj)
            obj.m_value = atan2(value(obj.y), value(obj.x));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('atan2(y, x)\n');
            
            begin_child(obj);
            draw(obj.y);
            end_child(obj);
            
            last_child(obj);
            draw(obj.x);
            end_child(obj);
        end
        
        function [topsort, visited, n_elem] = topological_sort(obj, topsort, visited, n_elem)
            % Topological sort of expression graph by a dfs.
            
            % Initialize if second and third args are empty
            if nargin == 1
                % topsort = {};
                visited = [];
                topsort = cell(1e4, 1);
                n_elem = 0;
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited, n_elem] = topological_sort(obj.y, topsort, visited, n_elem);
            [topsort, visited, n_elem] = topological_sort(obj.x, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end