classdef ast_timeinterval < yop.ast_unary_expression
    
    methods
        function obj = ast_timeinterval(t0, tf, expr)
            value = expr.m_value;
            ival = struct('bool', true, 't0', t0, 'tf', tf);
            obj@yop.ast_unary_expression(true, ival);
            obj.dim = expr.dim;
            obj.t0 = t0;
            obj.tf = tf;
            obj.expr = expr;
        end
        
        function value = evaluate(obj)            
            value = evaluate(obj.expr);
        end
        
        function v = forward(obj)
            obj.m_value = obj.expr.m_value;
            v = obj.m_value;
        end
        
        function ast(obj)
            fprintf(['[', num2str(obj.id), ']:', ...
                'interval(t0, tf, expr)\n']);
            
            begin_child(obj);
            ast(obj.t0);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.tf);
            end_child(obj);
            
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
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end