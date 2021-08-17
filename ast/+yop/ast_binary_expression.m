classdef ast_binary_expression < yop.ast_expression
    
    properties
        lhs
        rhs
    end
    
    methods
        function obj = ast_binary_expression(lhs, rhs)
            obj@yop.ast_expression();
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
        
        function draw(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            draw(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            draw(obj.rhs);
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
            [topsort, visited]=topological_sort(obj.lhs, topsort, visited);
            [topsort, visited]=topological_sort(obj.rhs, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
        
    end
end