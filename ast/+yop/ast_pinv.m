classdef ast_pinv < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_pinv(A)
            obj@yop.ast_expression();
            obj.A = A;
            obj.dim = size(pinv(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = pinv(evaluate(obj.A));
        end
        
        function v = forward(obj)
            obj.m_value = mod(value(obj.A));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('pinv(A)\n');
            last_child(obj);
            draw(obj.A);
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
            [topsort, visited, n_elem] = topological_sort(obj.A, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end