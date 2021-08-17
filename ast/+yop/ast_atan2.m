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
            [topsort, visited] = topological_sort(obj.y, topsort, visited);
            [topsort, visited] = topological_sort(obj.x, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end