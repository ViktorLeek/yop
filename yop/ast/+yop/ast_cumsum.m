classdef ast_cumsum < yop.ast_expression
    properties
        m_expr
        m_d
        m_nargs
    end
    methods
        function obj = ast_cumsum(expr, d)
            if nargin == 1
                val = cumsum(expr.m_value);
                num = cumsum(expr.m_numval);
                d = [];
            else
                val = cumsum(expr.m_value, d);
                num = cumsum(expr.m_numval, d);
                d = [];
            end
            sz = size(num);
            reducible = all(expr.m_reducible) & true(sz);
            t0 = max(expr.m_t0(:)) * ones(sz);
            tf = min(expr.m_tf(:)) * ones(sz);
            obj@yop.ast_expression( ...
                val      , ... value
                num      , ... numval
                t0       , ... t0
                tf       , ... tf
                false(sz), ... der
                reducible, ... reducible
                zeros(sz), ... type
                zeros(sz) ... typeid
                );
            obj.m_expr  = expr;
            obj.m_d     = d;
            obj.m_nargs = nargin;
        end
        
        function ast(obj)
            switch obj.m_nargs
                case 1
                    fprintf('cumsum(A)\n');
                    
                    last_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                case 2
                    fprintf('cumsum(A, dim)\n');
                    
                    begin_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_d);
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
            if ~isempty( find(visited == obj.m_id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.m_id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_expr, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_d, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end