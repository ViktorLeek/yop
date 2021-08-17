classdef ast_unary_expression < yop.ast_expression
    
    properties
        expr
    end
    
    methods
        function obj = ast_unary_expression(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
        end
        
        function draw(obj)
            fprintf([obj.name, '(expr)\n']);
            last_child(obj);
            draw(obj.expr);
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
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
        
    end
end