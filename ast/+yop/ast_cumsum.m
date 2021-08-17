classdef ast_cumsum < yop.ast_expression
    properties
        A
        d
        nargs
    end
    methods
        function obj = ast_cumsum(A, d)
            obj@yop.ast_expression();
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.A = A;
                    obj.dim = size(cumsum(ones(size(A))));
                    
                case 2
                    obj.A = A;
                    obj.d = d;
                    obj.dim = size(cumsum(ones(size(A)), d));
                    
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = cumsum(evaluate(obj.A));
                    
                case 2
                    value = cumsum(evaluate(obj.A), evaluate(obj.d));
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 1
                    obj.m_value = cumsum(value(obj.A));
                    
                case 2
                    obj.m_value = cumsum(value(obj.A), value(obj.d));
                    
            end
            v = obj.m_value;
        end
        
        function draw(obj)
            switch obj.nargs
                case 1
                    fprintf('cumsum(A)\n');
                    
                    last_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('cumsum(A, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.d);
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
            [topsort, visited, n_elem]=topological_sort(obj.A, topsort, visited, n_elem);
            [topsort, visited, n_elem]=topological_sort(obj.d, topsort, visited, n_elem);
            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
        
    end
end