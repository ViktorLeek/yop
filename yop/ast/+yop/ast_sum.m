classdef ast_sum < yop.ast_expression
    properties
        m_expr
        m_opt1
        m_opt2
        m_opt3
        m_nargs
    end
    methods
        function obj = ast_sum(expr, opt1, opt2, opt3)
            % ast_sum
            
            % This only works for cases where the optional arguments are
            % not yop.ast_node
            switch nargin
                case 1
                    val = sum(expr.m_value);
                    num = sum(expr.m_numval);
                    opt1=[]; opt2=[]; opt3=[];
                case 2
                    val = sum(expr.m_value, opt1);
                    num = sum(expr.m_numval, opt1);
                    opt2=[]; opt3=[];
                    
                case 3
                    val = sum(expr.m_value, opt1, opt2);
                    num = sum(expr.m_numval, opt1, opt2);
                    opt3=[];
                    
                case 4                    
                    val = sum(expr.m_value, opt1, opt2, opt3);
                    num = sum(expr.m_numval, opt1, opt2, opt3);
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
                false(sz), ... isder
                reducible, ... isreducible
                zeros(sz), ... type
                zeros(sz) ... typeid
                )
            obj.m_expr = expr;
            obj.m_opt1 = opt1;
            obj.m_opt2 = opt2;
            obj.m_opt3 = opt3;
            obj.m_nargs = nargin;
        end
        
        function ast(obj)
            
            switch obj.m_nargs
                case 1
                    fprintf('sum(expr)\n');
                    last_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                case 2
                    fprintf('sum(expr, opt1)\n');
                    
                    begin_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_opt1);
                    end_child(obj);
                    
                case 3
                    fprintf('sum(expr, opt1, opt2)\n');
                    
                    begin_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_opt1);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_opt2);
                    end_child(obj);
                    
                case 4
                    fprintf('sum(expr, opt1, opt2, opt3)\n');
                    
                    begin_child(obj);
                    ast(obj.m_expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_opt1);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.m_opt2);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.m_opt3);
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
                topological_sort(obj.m_expr, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_opt1, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_opt2, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_opt3, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end