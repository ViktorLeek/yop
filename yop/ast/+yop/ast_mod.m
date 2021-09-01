classdef ast_mod < yop.ast_expression
    properties
        a
        m
    end
    methods
        function obj = ast_mod(a, m)
            obj@yop.ast_expression();
            obj.a = a;
            obj.m = m;
            obj.dim = size(mod(ones(size(a)), ones(size(m))));
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            if all(isa_numeric(obj.a)) && all(isa_numeric(obj.m))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function boolv = is_transcription_invariant(obj)
            if all(is_transcription_invariant(obj.a)) && ...
                    all(is_transcription_invariant(obj.m))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function obj = set_pred(obj)
            add_pred(obj.a, obj);
            add_pred(obj.m, obj);
        end
        
        function value = evaluate(obj)
            value = mod(evaluate(obj.a), evaluate(obj.m));
        end
        
        function v = forward(obj)
            obj.m_value = mod(value(obj.a), value(obj.m));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('mod(a, m)\n');
            
            begin_child(obj);
            draw(obj.a);
            end_child(obj);
            
            last_child(obj);
            draw(obj.m);
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
                topological_sort(obj.a, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m, visited, topsort, n_elem);

            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end