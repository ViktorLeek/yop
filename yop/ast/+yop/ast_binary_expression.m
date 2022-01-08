classdef ast_binary_expression < yop.ast_expression
    
    properties
        m_lhs
        m_rhs
    end
    
    methods
        function obj = ast_binary_expression(value, numval, t0, tf, isder, isreducible, type, typeid, lhs, rhs)
            obj@yop.ast_expression(value, numval, t0, tf, isder, isreducible, type, typeid)
            obj.m_lhs = lhs;
            obj.m_rhs = rhs;
        end
        
        function ast(obj)
            fprintf([obj.m_name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            ast(obj.m_lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m_rhs);
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
                topological_sort(obj.m_lhs, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.m_rhs, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end