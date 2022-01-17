classdef ast_linspace < yop.ast_expression
    properties
        m_x1
        m_x2
        m_n
    end
    methods
        function obj = ast_linspace(x1, x2, n)
            val = linspace(value(x1), value(x2), n);
            num = linspace(numval(x1), numval(x2), n);
            sz  = size(num);
            tmp = isa_reducible(x1) & isa_reducible(x2);
            reducible = all(tmp(:)) & true(sz);
            t0_1 = get_t0(x1);
            t0_2 = get_t0(x2);
            tf_1 = get_tf(x1);
            tf_2 = get_tf(x2);
            t0 = max([t0_1(:); t0_2(:)]) * ones(sz);
            tf = min([tf_1(:); tf_2(:)]) * ones(sz);
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
            
            obj.m_x1 = x1;
            obj.m_x2 = x2;
            obj.m_n = n;
            obj.dim = [1, n];
        end
        
        function ast(obj)
            fprintf('linspace(x1, x2, n)\n');
            
            begin_child(obj);
            ast(obj.m_x1);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.m_x2);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_n);
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
                topological_sort(obj.m_x1, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_x2, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_n, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end