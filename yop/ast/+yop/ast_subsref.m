classdef ast_subsref < yop.ast_expression
    
    properties
        m_expr
        m_s
    end
    
    methods
        function obj = ast_subsref(expr, s)
            % Casadi is a bit picky, so logical indices must be converted
            % to indices.
            for k=1:length(s.subs)
                if islogical(s.subs{k})
                    s.subs{k} = find(s.subs{k});
                end
            end
            
            obj@yop.ast_expression( ...
                subsref(expr.m_value    , s), ...
                subsref(expr.m_numval   , s), ...
                subsref(expr.m_t0       , s), ...
                subsref(expr.m_tf       , s), ...
                subsref(expr.m_der      , s), ...
                subsref(expr.m_reducible, s), ...
                subsref(expr.m_type     , s), ...
                subsref(expr.m_typeid   , s) ...
                );
            obj.m_expr = expr;
            obj.m_s = s;
        end
        
        function ast(obj)
            fprintf('subsref(node, s)\n');
            
            begin_child(obj);
            ast(obj.m_expr);
            end_child(obj);
            
            str = [];
            for k=1:length(obj.m_s.subs)
                idx = obj.m_s.subs{k};
                if size(idx,1) > size(idx,2)
                    idx = idx';
                end
                str = [str, '[', num2str(idx), '], '];
            end
            last_child(obj);
            fprintf(['{', str(1:end-2), '}\n']);
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
            if ~isempty( find(visited == obj.m_id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.m_id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_expr, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_s, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end