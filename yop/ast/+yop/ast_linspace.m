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
        
        function boolv = isa_numeric(obj)
            if all(isa_numeric(obj.x1)) && all(isa_numeric(obj.x2)) ...
                    && all(isa_numeric(obj.n))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.x1)) ...
                    && all(isa_reducible(obj.x2)) ...
                    && all(isa_reducible(obj.n))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
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
                topological_sort(obj.x1, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.x2, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.n, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end