classdef ast_cross < yop.ast_expression
    properties
        A
        B
        d
        nargs
    end
    methods
        function obj = ast_cross(A, B, d)
            obj@yop.ast_expression();
            obj.nargs = nargin;
            switch nargin
                case 2
                    obj.A = A;
                    obj.B = B;
                    obj.dim = size(cross(ones(size(A)), ones(size(B))));
                    
                case 3
                    obj.A = A;
                    obj.B = B;
                    obj.d = d;
                    obj.dim = size(cross(ones(size(A)), ones(size(B)), d));
                    
                otherwise
                    error('Wrong number of arguments provided')
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 2
                    value = cross(evaluate(obj.A), evaluate(obj.B));
                case 3
                    value = cross(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 2
                    obj.m_value = cross(value(obj.A), value(obj.B));
                    
                case 3
                    obj.m_value = cross(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d) ...
                        );
            end            
            v = obj.m_value;
        end
        
        function draw(obj)
            switch obj.nargs
                case 2
                    fprintf('cross(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('cross(A, B, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
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
            [topsort, visited, n_elem] = topological_sort(...
                obj.A, ...
                topsort, ...
                visited, ...
                n_elem ...
                );
            
            [topsort, visited, n_elem] = topological_sort(...
                obj.B, ...
                topsort, ...
                visited, ...
                n_elem ...
                );
            
            % If d is empty, this call should be dispatched to function
            % topological sort, as opposed to the method, so it should
            % simply return the topsort and visited.
            [topsort, visited, n_elem] = topological_sort(...
                obj.d, ...
                topsort, ...
                visited, ...
                n_elem ...
                );
            

            
            % append self to sort
            n_elem = n_elem + 1;
            topsort{n_elem} = obj;

        end
    end
end