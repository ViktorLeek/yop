classdef ast_dot < yop.ast_expression
    properties
        A
        B
        d
        nargs
    end
    methods
        function obj = ast_dot(A, B, d)
            obj@yop.ast_expression();
            obj.nargs = nargin;
            switch nargin
                case 2
                    obj.A = A;
                    obj.B = B;
                    obj.dim = size(dot(ones(size(A)), ones(size(B))));
                    
                case 3
                    obj.A = A;
                    obj.B = B;
                    obj.d = d;
                    obj.dim = size(dot(ones(size(A)), ones(size(B)), d));
                    
                otherwise
                    error('Wrong number of arguments provided')
            end
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 2
                    value = dot(evaluate(obj.A), evaluate(obj.B));
                    
                case 3
                    value = dot(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 2
                    obj.m_value = dot(value(obj.A), value(obj.B));
                    
                case 3
                    obj.m_value = dot(...
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
                    fprintf('dot(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('dot(A, B, dim)\n');
                    
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
        
        function [topsort, visited] = topological_sort(obj, topsort, visited)
            % Topological sort of expression graph by a dfs.
            
            % Initialize if second and third args are empty
            if nargin == 1
                topsort = {};
                visited = [];
            end
            
            % only visit every node once
            if ~isempty( find(visited == obj.id, 1) )
                return;
            end
            
            % Mark node as visited
            visited = [visited, obj.id];
            
            % Visit child
            [topsort, visited] = topological_sort(obj.A, topsort, visited);
            [topsort, visited] = topological_sort(obj.B, topsort, visited);
            [topsort, visited] = topological_sort(obj.d, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end