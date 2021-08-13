classdef ast_cross < yop.ast_node
    properties
        A
        B
        d
        nargs
    end
    methods
        function obj = ast_cross(A, B, d)
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
        
        function ast(obj)
            switch obj.nargs
                case 2
                    fprintf('cross(A, B)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('cross(A, B, dim)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.d);
                    end_child(obj);
                    
            end
            
        end
    end
end