classdef ast_min < yop.ast_expression
    properties
        A
        B
        d
        flag
        nargs
    end
    methods
        function obj = ast_min(A, B, d, flag)
            obj@yop.ast_expression();
            obj.A = A;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(min(ones(size(A))));
                    
                case 2
                    obj.B = B;
                    if isa(B, 'yop.node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp));
                    
                case 3
                    obj.B = B;
                    obj.d = d;
                    if isa(B, 'yop.node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d));
                    
                case 4
                    obj.B = B;
                    obj.d = d;
                    obj.flag = flag;
                    if isa(B, 'yop.node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d, flag));     
            end
            
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = min(evaluate(obj.A));
                    
                case 2
                    value = min(evaluate(obj.A), evaluate(obj.B));
                    
                case 3
                    value = min(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
                    
                case 4
                    value = min(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d), ...
                        evaluate(obj.flag) ...
                        );
            end
        end
        
        function v = forward(obj)
            switch obj.nargs
                case 1
                    obj.m_value = min(value(obj.A));
                    
                case 2
                    obj.m_value = min(value(obj.A), value(obj.B));
                    
                case 3
                    obj.m_value = min(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d) ...
                        );
                    
                case 4
                    obj.m_value = min(...
                        value(obj.A), ...
                        value(obj.B), ...
                        value(obj.d), ...
                        value(obj.flag) ...
                        );
            end         
            v = obj.m_value;
        end
        
        function draw(obj)
            switch obj.nargs
                case 1
                    fprintf('min(A)\n');
                    last_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('min(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('min(A, B, dim)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.d);
                    end_child(obj);
                    
                case 4
                    fprintf('min(A, B, dim, flag)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                    begin_child(obj);
                    draw(obj.d);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.flag);
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
            [topsort, visited]=topological_sort(obj.A, topsort, visited);
            [topsort, visited]=topological_sort(obj.b, topsort, visited);
            [topsort, visited]=topological_sort(obj.d, topsort, visited);
            [topsort, visited]=topological_sort(obj.flag, topsort, visited);
            
            % append self to sort
            topsort = [topsort(:)', {obj}];
        end
    end
end