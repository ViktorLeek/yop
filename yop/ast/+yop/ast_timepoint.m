classdef ast_timepoint < yop.ast_expression
    properties
        m_timepoint
        m_expr
    end
    methods
        function obj = ast_timepoint(timepoint, expr)
            sz = size(expr);
            obj@yop.ast_expression( ...
                yop.cx('tp', sz(1), sz(2)), ... value
                expr.m_numval             , ... numval
                timepoint                 , ... t0
                timepoint                 , ... tf
                false(sz)                 , ... der
                true(sz)                  , ... reducible
                expr.m_type               , ... type
                expr.m_typeid              ... typeid
                );
            obj.m_timepoint = timepoint;
            obj.m_expr = expr;
        end
        
        function ast(obj)
            fprintf(['[', num2str(obj.m_id), ']:', ...
                'timepoint(timepoint, expr)\n']);
            
            begin_child(obj);
            ast(obj.m_timepoint);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_expr);
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
                topological_sort(obj.m_timepoint, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end