classdef ast_der < yop.ast_expression
    properties
        var
    end
    methods
        function obj = ast_der(var)
            obj@yop.ast_expression();
            obj.var = var;
            obj.dim = size(var);
        end
        
        % Overload all illegal operations here!! which should be most!
        
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
        
        function bool = is_differential(obj)
            bool = true;
        end
        
        function draw(obj)
            fprintf('der(var)\n');
            last_child(obj);
            draw(obj.var);
            end_child(obj);
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
                topological_sort(obj.var, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end