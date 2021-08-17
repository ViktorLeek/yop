classdef ast_alg < yop.ast_expression
    properties
        var
    end
    methods
        function obj = ast_alg(var)
            obj@yop.ast_expression();
            obj.var = var;
            obj.dim = size(var);
        end
        
        function value = evaluate(obj)
            value = evaluate(obj.var);
        end
        
        function value = forward(obj)
            error([...
                '[yop] Error: Cannot forward evaluate an ', ...
                class(obj), ...
                ' node.' ...
                ]);
        end
        
        function draw(obj)
            fprintf('alg(var)\n');
            last_child(obj);
            draw(obj.var);
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
            [topsort, visited]=topological_sort(obj.var, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end