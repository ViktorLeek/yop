classdef ast_int < yop.ast_expression
    properties
        expr
    end
    methods
        function obj = ast_int(expr)
            isival = false; % Always finite integrals (for now at least)
            obj@yop.ast_expression(isival);
            obj.expr = expr;
            obj.dim = size(expr);
            obj.m_value = yop.cx('int', size(expr,1), size(expr,2));
        end
        
        function val = numval(obj)
            val = numval(obj.expr);
        end
        
        function boolv = isa_reducible(obj)
            boolv = true(size(obj));
        end
        
        function value = evaluate(obj)
            % Evaluates like an ast_variable. The reason is that it does
            % not have the semantics of a timevarying expression, so it
            % needs to be parameterized separately. Before parametrization
            % has occured, the idea is that an MX/sym variable (of the
            % right size) is used to represent its value.
            value = obj.m_value;
        end
        
        function v = forward(obj)            
            % Evaluates like an ast_variable. The reason is that it does
            % not have the semantics of a timevarying expression, so it
            % needs to be parameterized separately. Before parametrization
            % has occured, the idea is that an MX/sym variable (of the
            % right size) is used to represent its value.
            v = obj.m_value;
        end
        
        function ast(obj)
            fprintf('int(expr)\n');
            last_child(obj);
            ast(obj.expr);
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
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;
        end 
    end
end