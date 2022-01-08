classdef ast_cross < yop.ast_expression
    properties
        m_A
        m_B
        m_d
        m_nargs
    end
    methods
        function obj = ast_cross(A, B, d)
            if nargin == 2
                val = cross(value(A), value(B));
                num = cross(numval(A), numval(B));
                d = [];
            else
                val = cross(value(A), value(B), d);
                num = cross(numval(A), numval(B), d);
            end
            
            sz  = size(num);
            reducible = all(isa_reducible(A) & isa_reducible(B)) & true(sz);
            t0_A = get_t0(A);
            t0_B = get_t0(B);
            tf_A = get_tf(A);
            tf_B = get_tf(B);
            t0 = max([t0_A(:); t0_B(:)]) * ones(sz);
            tf = min([tf_A(:); tf_B(:)]) * ones(sz);
            
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
            obj.m_A = A;
            obj.m_B = B;
            obj.m_d = d;
            obj.m_nargs = nargin;
        end
        
        function ast(obj)
            if obj.m_nargs == 2
                fprintf('cross(A, B)\n');
                
                begin_child(obj);
                ast(obj.m_A);
                end_child(obj);
                
                last_child(obj);
                ast(obj.m_B);
                end_child(obj);
                
            else
                fprintf('cross(A, B, dim)\n');
                
                begin_child(obj);
                ast(obj.m_A);
                end_child(obj);
                
                begin_child(obj);
                ast(obj.m_B);
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
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_A, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_B, visited, topsort, n_elem);
            
            % If d is empty, this call should be dispatched to function
            % topological_sort, as opposed to the method, so it should
            % simply return the topsort and visited.
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_d, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end