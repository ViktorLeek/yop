classdef ast_interp1 < yop.ast_expression
    properties
        m_x
        m_v
        m_xq
    end
    methods
        function obj = ast_interp1(x, v, xq)
            lut = casadi.interpolant('interp1', 'bspline', {x}, v);
            val = lut(xq.m_value);
            num = full(lut(x(1)*ones(size(xq))));
            sz = size(num);
            reducible = false(sz);
            t0 = max(xq.m_t0(:)) * ones(sz);
            tf = min(xq.m_tf(:)) * ones(sz);
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
            obj.m_x  = x;
            obj.m_v  = v;
            obj.m_xq = xq;
        end
        
        function ast(obj)
            
            begin_child(obj);
            ast(obj.m_x);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.m_v);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_xq);
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
                topological_sort(obj.m_xq, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end