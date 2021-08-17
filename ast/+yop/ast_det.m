classdef ast_det < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_det(A)
            obj@yop.ast_expression();
            obj.A = A;
            obj.dim = size(det(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = det(evaluate(obj.A));
        end
        
        function v = forward(obj)
            obj.m_value = det(value(obj.A));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('det(A)\n');
            last_child(obj);
            draw(obj.A);
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
            [topsort, visited]=topological_sort(obj.A, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end