classdef ast_int < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_int(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(expr);
        end
        
        function value = evaluate(obj)
            error('[Yop] Error: Unable to evaluate node.');
        end
        
        function v = forward(obj)
            error('[Yop] Error: Unable to forward evaluate node.');
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
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
        function draw(obj)
            fprintf('int(expr)\n');
            last_child(obj);
            draw(obj.expr);
            end_child(obj);
        end
    end
end