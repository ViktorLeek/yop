classdef ast_dot < yop.ast_node
    properties
        A
        B
        d
        nargs
    end
    methods
        function obj = ast_dot(A, B, d)
            switch nargin
                case 2
                    obj.nargs = 2;
                    obj.A = A;
                    obj.B = B;
                    obj.dim = size(dot(ones(size(A)), ones(size(B))));
                    
                case 3
                    obj.nargs = 3;
                    obj.A = A;
                    obj.B = B;
                    obj.d = d;
                    obj.dim = size(dot(ones(size(A)), ones(size(B)), d));
                    
                otherwise
                    error('Wrong number of arguments provided')
            end
        end
        function ast(obj)
            switch obj.nargs
                case 2
                    fprintf('dot(A, B)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('dot(A, B, dim)\n');
                    
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