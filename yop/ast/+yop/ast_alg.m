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
        
        % Overload all illegal operations here!! which should be most!
        % only allow 0 = alg(expr)!
        
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
        
        function obj = set_pred(obj)
            add_pred(obj.var, obj);
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            error('[Yop] Error: Method not implemented.');
        end
        
        function draw(obj)
            fprintf('alg(var)\n');
            last_child(obj);
            draw(obj.var);
            end_child(obj);
        end
        
        function [topsort, n_elem, visited] = ...
                topological_sort(obj, visited, topsort, n_elem)
            % Topological sort of expression graph by a dfs.
            
            switch nargin
                case 1
                    % Start new sort: topological_sort(obj)
                    visited = [];
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                case 2
                    % Semi-warm start: topological_sort(obj, visited)
                    % In semi-warm start 'visited' is already provided, but 
                    % no elements are sorted. This is for instance useful 
                    % for finding all variables in a number of expressions 
                    % that are suspected to contain common subexpressions.
                    topsort = ...
                        cell(yop.constants().topsort_preallocation_size, 1);
                    n_elem = 0;
                    
                otherwise
                    % Pass
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end
