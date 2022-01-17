classdef ast_max < yop.ast_expression
    properties
        m_A
        m_B
        m_d
        m_flag
        m_nargs
    end
    methods
        function obj = ast_max(A, B, d, flag)
            switch nargin
                case 1
                    val = max(A.m_value);
                    num = max(A.m_numval);
                    sz  = size(num);
                    red = all(A.m_reducible(:)) * ones(sz);
                    t0  = A.m_t0;
                    tf  = A.m_tf;
                    B=[]; d=[]; flag=[];
                    
                case 2
                    val = max(value(A), value(B));
                    num = max(numval(A), numval(B));
                    sz  = size(num);
                    tmp = isa_reducible(A) & isa_reducible(B);
                    red = all(tmp(:)) & true(sz);
                    t0_A = get_t0(A);
                    t0_B = get_t0(B);
                    tf_A = get_tf(A);
                    tf_B = get_tf(B);
                    t0 = max([t0_A(:); t0_B(:)]) * ones(sz);
                    tf = min([tf_A(:); tf_B(:)]) * ones(sz);
                    d=[]; flag=[];
                    
                case 3
                    val = max(value(A), value(B), d);
                    num = max(numval(A), numval(B), d);
                    sz  = size(num);
                    tmp = isa_reducible(A) & isa_reducible(B);
                    red = all(tmp(:)) & true(sz);
                    t0_A = get_t0(A);
                    t0_B = get_t0(B);
                    tf_A = get_tf(A);
                    tf_B = get_tf(B);
                    t0 = max([t0_A(:); t0_B(:)]) * ones(sz);
                    tf = min([tf_A(:); tf_B(:)]) * ones(sz);
                    flag=[];
                    
                case 4
                    val = max(value(A), value(B), d, flag);
                    num = max(numval(A), numval(B), d, flag);
                    sz  = size(num);
                    tmp = isa_reducible(A) & isa_reducible(B);
                    red = all(tmp(:)) & true(sz);
                    t0_A = get_t0(A);
                    t0_B = get_t0(B);
                    tf_A = get_tf(A);
                    tf_B = get_tf(B);
                    t0 = max([t0_A(:); t0_B(:)]) * ones(sz);
                    tf = min([tf_A(:); tf_B(:)]) * ones(sz);
            end
            
            obj@yop.ast_expression( ...
                val      , ... value
                num      , ... numval
                t0       , ... t0
                tf       , ... tf
                false(sz), ... der
                red      , ... reducible
                zeros(sz), ... type
                zeros(sz) ... typeid
                );
            obj.m_A = A;
            obj.m_B = B;
            obj.m_d = d;
            obj.m_flag  = flag;
            obj.m_nargs = nargin;
        end

        
        function ast(obj)
            switch obj.m_nargs
                case 1
                    fprintf('max(A)\n');
                    last_child(obj);
                    ast(obj.m_A);
                    end_child(obj);
                    
                case 2
                    fprintf('max(A, B)\n');
                    
                    begin_child(obj);
                    ast(obj.m_A);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_B);
                    end_child(obj);
                    
                case 3
                    fprintf('max(A, B, dim)\n');
                    
                    begin_child(obj);
                    ast(obj.m_A);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_B);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_d);
                    end_child(obj);
                    
                case 4
                    fprintf('max(A, B, dim, flag)\n');
                    
                    begin_child(obj);
                    ast(obj.m_A);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_B);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_d);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_flag);
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
                topological_sort(obj.m_A, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_B, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_d, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_flag, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end