classdef ast_linspace < yop.ast_expression
    properties
        x1
        x2
        n
    end
    methods
        function obj = ast_linspace(x1, x2, n)
            obj@yop.ast_expression();
            obj.x1 = x1;
            obj.x2 = x2;
            obj.n = n;
            obj.dim = [1, n];
        end
        
        function value = evaluate(obj)
            value = linspace(...
                evaluate(obj.x1), ...
                evaluate(obj.x2), ...
                evaluate(obj.n) ...
                );
        end
        
        function v = forward(obj)
            obj.m_value = linspace(... 
                value(obj.x1), ...
                value(obj.x2), ...
                value(obj.n) ...
                );
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('linspace(x1, x2, n)\n');
            
            begin_child(obj);
            draw(obj.x1);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.x2);
            end_child(obj);
            
            last_child(obj);
            draw(obj.n);
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
            [topsort, visited] = topological_sort(obj.x1, topsort, visited);
            [topsort, visited] = topological_sort(obj.x2, topsort, visited);
            [topsort, visited] = topological_sort(obj.n, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end