classdef ast_binary_expression < yop.ast_expression
    
    properties
        lhs
        rhs
    end
    
    methods
        function obj = ast_binary_expression(lhs, rhs)
            obj@yop.ast_expression(is_ival(lhs) || is_ival(rhs));
            obj.lhs = lhs;
            obj.rhs = rhs;
        end
        
        function val = numval(obj)
            val_lhs = numval(obj.lhs);
            val_rhs = numval(obj.rhs);
            switch obj.name
                case 'plus'
                    val = plus(val_lhs, val_rhs);
                case 'minus'
                    val = minus(val_lhs, val_rhs);
                case 'mtimes'
                    val = mtimes(val_lhs, val_rhs);
                case 'mpower'
                    val = mpower(val_lhs, val_rhs);
                case 'mrdivide'
                    val = mrdivide(val_lhs, val_rhs);
                case 'power'
                    val = power(val_lhs, val_rhs);
                case 'times'
                    val = times(val_lhs, val_rhs);
                case 'rdivide'
                    val = rdivide(val_lhs, val_rhs);
                case 'mldivide'
                    val = mldivide(val_lhs, val_rhs);
                case 'ldivide'
                    val = ldivide(val_lhs, val_rhs);
            end            
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.lhs)) && ...
                    all(isa_reducible(obj.rhs))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function ast(obj)
            fprintf([obj.name, '(lhs, rhs)\n']);
            
            begin_child(obj);
            ast(obj.lhs);
            end_child(obj);
            
            last_child(obj);
            ast(obj.rhs);
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
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, n_elem, visited] = ...
                topological_sort(obj.lhs, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.rhs, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end