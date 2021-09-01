classdef ast_inv < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_inv(A)
            obj@yop.ast_expression();
            obj.A = A;
            obj.dim = size(inv(ones(size(A))));
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            boolv = isa_numeric(obj.A);
        end
        
        function boolv = is_transcription_invariant(obj)
            boolv = is_transcription_invariant(obj.expr);
        end
        
        function obj = set_pred(obj)
            add_pred(obj.A, obj);
        end
        
        function value = evaluate(obj)
            value = inv(evaluate(obj.A));
        end
        
        function v = forward(obj)
            obj.m_value = inv(value(obj.A));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('inv(A)\n');
            last_child(obj);
            draw(obj.A);
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
                topological_sort(obj.A, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end