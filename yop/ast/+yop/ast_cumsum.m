classdef ast_cumsum < yop.ast_expression
    properties
        A
        d
        nargs
    end
    methods
        function obj = ast_cumsum(A, d)
            obj@yop.ast_expression();
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.A = A;
                    obj.dim = size(cumsum(ones(size(A))));
                    
                case 2
                    obj.A = A;
                    obj.d = d;
                    obj.dim = size(cumsum(ones(size(A)), d));
                    
            end
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            boolv = isa_numeric(obj.expr);
            switch obj.nargs
                case 1
                    if all(isa_numeric(obj.A))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
                case 2
                    if all(isa_numeric(obj.A)) && all(isa_numeric(obj.d))
                        boolv = true(size(obj));
                    else
                        boolv = false(size(obj));
                    end
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = cumsum(evaluate(obj.A));
                    
                case 2
                    value = cumsum(evaluate(obj.A), evaluate(obj.d));
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 1
                    obj.m_value = cumsum(value(obj.A));
                    
                case 2
                    obj.m_value = cumsum(value(obj.A), value(obj.d));
                    
            end
            v = obj.m_value;
        end
        
        function draw(obj)
            switch obj.nargs
                case 1
                    fprintf('cumsum(A)\n');
                    
                    last_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('cumsum(A, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.d);
                    end_child(obj);
            end
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
                topological_sort(obj.A, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.d, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end