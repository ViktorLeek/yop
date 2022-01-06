classdef ast_atan2 < yop.ast_expression
    properties
        y
        x
    end
    methods
        function obj = ast_atan2(y, x)
            obj@yop.ast_expression(is_ival(x) || is_ival(y));
            obj.y = y;
            obj.x = x;
            obj.dim = size(atan2(ones(size(y)), ones(size(x))));
        end
        
        function boolv = isa_numeric(obj)
            if all(isa_numeric(obj.y) & isa_numeric(obj.x))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function boolv = is_transcription_invariant(obj)
            if all(is_transcription_invariant(obj.y)) && ...
                    all(is_transcription_invariant(obj.x))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function value = evaluate(obj)
            value = atan2(evaluate(obj.y), evaluate(obj.x));
        end
        
        function v = forward(obj)
            obj.m_value = atan2(value(obj.y), value(obj.x));
            v = obj.m_value;
        end
        
        function draw(obj)
            fprintf('atan2(y, x)\n');
            
            begin_child(obj);
            draw(obj.y);
            end_child(obj);
            
            last_child(obj);
            draw(obj.x);
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
                topological_sort(obj.y, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.x, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end