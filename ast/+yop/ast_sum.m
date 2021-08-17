classdef ast_sum < yop.ast_expression
    properties
        expr
        opt1
        opt2
        opt3
        nargs
    end
    methods
        function obj = ast_sum(expr, opt1, opt2, opt3)
            % ast_sum
            
            % This only works for cases where the optional arguments are
            % not yop.node
            obj@yop.ast_expression();
            obj.expr = expr;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(sum(ones(size(expr))));
                    
                case 2
                    obj.opt1 = opt1;
                    obj.dim = size(sum(ones(size(expr)), opt1));
                    
                case 3
                    obj.opt1 = opt1;
                    obj.opt2 = opt2;
                    obj.dim = size(sum(ones(size(expr)), opt1, opt2));
                    
                case 4                    
                    obj.opt1 = opt1;
                    obj.opt2 = opt2;
                    obj.opt3 = opt3;
                    obj.dim = size(sum(ones(size(expr)), opt1, opt2, opt3));
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = sum(evaluate(obj.expr));
                    
                case 2
                    value = sum(evaluate(obj.expr), obj.opt1);
                    
                case 3
                    value = sum(evaluate(obj.expr), obj.opt1, obj.opt2);
                    
                case 4
                    value = sum(...
                        evaluate(obj.expr), ...
                        obj.opt1, ...
                        obj.opt2, ...
                        obj.opt3 ...
                        );
                    
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 1
                    obj.m_value = sum(value(obj.expr));
                    
                case 2
                    obj.m_value = sum(value(obj.expr), obj.opt1);
                    
                case 3
                    obj.m_value = sum(value(obj.expr), obj.opt1, obj.opt2);
                    
                case 4
                    obj.m_value = sum( ...
                        value(obj.expr), ...
                        obj.opt1, ...
                        obj.opt2, ...
                        obj.opt3 ...
                        ); 
            end
            v = obj.m_value;
        end
        
        function draw(obj)
            
            switch obj.nargs
                case 1
                    fprintf('sum(expr)\n');
                    last_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                case 2
                    fprintf('sum(expr, opt1)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                case 3
                    fprintf('sum(expr, opt1, opt2)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt2);
                    end_child(obj);
                    
                case 4
                    fprintf('sum(expr, opt1, opt2, opt3)\n');
                    
                    begin_child(obj);
                    draw(obj.expr);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt1);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.opt2);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.opt3);
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
            [topsort, visited, n_elem]=topological_sort(obj.expr, topsort, visited, n_elem);
            [topsort, visited, n_elem]=topological_sort(obj.opt1, topsort, visited, n_elem);
            [topsort, visited, n_elem]=topological_sort(obj.opt2, topsort, visited, n_elem);
            [topsort, visited, n_elem]=topological_sort(obj.opt3, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end