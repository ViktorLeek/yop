classdef ast_norm < yop.ast_expression
    properties
        expr
        p
        nargs
    end
    methods
        function obj = ast_norm(expr, p)
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.nargs = nargin;
            obj.dim = [1, 1];
            if nargin==2
                obj.p = p;
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
        
        function [topsort, visited, n_elem] = topological_sort(obj, topsort, visited, n_elem)
            % Topological sort of expression graph by a dfs.
            
            % Initialize if second and third args are empty
            if nargin == 1
                % topsort = {};
                visited = [];
                topsort = cell(1e4, 1);
                n_elem = 0;
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited, n_elem] = topological_sort(obj.expr, topsort, visited, n_elem);
            [topsort, visited, n_elem] = topological_sort(obj.p, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end