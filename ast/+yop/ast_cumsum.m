classdef ast_cumsum < yop.ast_node
    properties
        A
        d
        nargs
    end
    methods
        function obj = ast_dot(A, d)
            switch nargin
                case 1
                    obj.nargs = 1;
                    obj.A = A;
                    obj.dim = size(cumsum(ones(size(A))));
                    
                case 2
                    obj.nargs = 2;
                    obj.A = A;
                    obj.d = d;
                    obj.dim = size(cumsum(ones(size(A)), d));
                    
            end
        end
        function ast(obj)
            switch obj.nargs
                case 1
                    fprintf('cumsum(A, B)\n');
                    
                    last_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('cumsum(A, dim)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.d);
                    end_child(obj);
                    
            end
            
        end
    end
end