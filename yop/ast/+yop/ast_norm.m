classdef ast_norm < yop.ast_expression
    properties
        expr
        p
        nargs
    end
    methods
        function obj = ast_norm(expr, p)
            obj@yop.ast_expression(is_ival(expr));
            obj.expr = expr;
            obj.nargs = nargin;
            obj.dim = [1, 1];
            if nargin==2
                obj.p = p;
            end
        end
        
        function val = numval(obj)
            if obj.nargs == 1
                val = norm(numval(obj.expr));
            else
                val = norm(numval(obj.expr), numval(obj.p));
            end
        end
        
        function boolv = isa_reducible(obj)
            if all(isa_reducible(obj.expr))
                boolv = true(size(obj));
            else
                boolv = false(size(obj));
            end
        end
        
        function value = evaluate(obj)
            if obj.nargs == 1
                value = norm(evaluate(obj.expr));
            else
                value = norm(evaluate(obj.expr), evaluate(obj.p));
            end
        end
        
        function v = forward(obj)
            if obj.nargs == 1
                obj.m_value = norm(value(obj.expr));
            else
                obj.m_value = norm(value(obj.expr), value(obj.p));
            end
            v = obj.m_value;
        end
        
        function draw(obj)
            if obj.nargs == 1
                fprintf('norm(expr)\n');
                last_child(obj);
                draw(obj.expr);
                end_child(obj);
            else
                fprintf('norm(expr, p)\n');
                begin_child(obj);
                draw(obj.expr);
                end_child(obj);
                last_child(obj);
                draw(obj.p);
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
                topological_sort(obj.expr, visited, topsort, n_elem);
            
            [topsort, n_elem, visited] = ...
                topological_sort(obj.p, visited, topsort, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end