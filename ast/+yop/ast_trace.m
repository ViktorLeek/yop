classdef ast_trace < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_trace(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(trace(ones(size(expr))));
        end
        
        function value = evaluate(obj)
            value = trace(evaluate(obj.expr));
        end
        
        function v = forward(obj)
            obj.m_value = trace(value(obj.expr));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('trace(expr)\n');
            last_child(obj);
            draw(obj.expr);
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
            [topsort, visited, n_elem]=topological_sort(obj.expr, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end