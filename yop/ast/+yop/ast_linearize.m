classdef ast_linearize < yop.ast_expression
    
    properties
        m_f
        m_x
        m_x0
    end
    
    methods
        function obj = ast_linearize(f, x, x0)
            val = linearize(value(f), value(x), value(x0));
            sz = size(val);
            num = zeros(sz);
            obj@yop.ast_expression( ...
                val      , ... value
                num      , ... numval
                -inf(sz) , ... t0
                +inf(sz) , ... tf
                false(sz), ... der
                false(sz), ... reducible
                zeros(sz), ... type
                zeros(sz) ... typeid
                );
            obj.m_f  = f;
            obj.m_x  = x;
            obj.m_x0 = x0;
        end
        
        function ast(obj)
            fprintf('linearize(f, x, x0)\n');
            
            begin_child(obj);
            ast(obj.m_f);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.m_x);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_x0);
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
                topological_sort(obj.m_f, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_x, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_x0, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
    
end