classdef ast_transpose < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_transpose(expr)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.dim = size(transpose(ones(size(expr))));
        end
        
        function boolv = isa_numeric(obj)
            % Potentially very slow. If it turns out to be too slow an
            % alternative solution, such as a DFS can be used.
            boolv = transpose(isa_numeric(obj.expr));
        end
        
        function [bool, tp] = isa_timepoint(obj)
            [bool, tp] = isa_timepoint(obj.expr);
            bool = transpose(bool);
            tp = transpose(tp);
        end
        
        function boolv = isa_der(obj)
            boolv = transpose(isa_der(obj.expr));
        end
        
        function [bool, id] = isa_variable(obj)
            [bool, id] = isa_variable(obj.expr);
            bool = transpose(bool);
            id = transpose(id);
        end
        
        
        function boolv = is_transcription_invariant(obj)
            boolv = transpose(is_transcription_invariant(obj.expr));
        end
        
        function obj = set_pred(obj)
            add_pred(obj.expr, obj);
        end
        
        function value = evaluate(obj)
            value = transpose(evaluate(obj.expr));
        end
        
        function v = forward(obj)
            obj.m_value = transpose(value(obj.expr));
            v = obj.m_value;
        end
    
        function draw(obj)
            fprintf('transpose(obj)\n');
            last_child(obj);
            draw(obj.expr);
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