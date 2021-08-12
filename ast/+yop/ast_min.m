classdef ast_min < yop.ast_node
    properties
        A
        B
        d
        flag
        n_args
    end
    methods
        function obj = ast_min(A, B, d, flag)
            obj.A = A;
            
            switch nargin
                case 1
                    obj.n_args = 1;
                    obj.dim = size(min(ones(size(A))));
                    
                case 2
                    obj.n_args = 2;
                    obj.B = B;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp));
                    
                case 3
                    obj.n_args = 3;
                    obj.B = B;
                    obj.d = d;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d));
                    
                case 4
                    obj.n_args = 4;
                    obj.B = B;
                    obj.d = d;
                    obj.flag = flag;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d, flag));     
            end
            
        end
        function ast(obj)
            switch obj.n_args
                case 1
                    fprintf('min(A)\n');
                    last_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('min(A, B)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('min(A, B, dim)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.d);
                    end_child(obj);
                    
                case 4
                    fprintf('min(A, B, dim, flag)\n');
                    
                    begin_child(obj);
                    ast(obj.A);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.B);
                    end_child(obj);
                    
                    begin_child(obj);
                    ast(obj.d);
                    end_child(obj);
                    
                    last_child(obj);
                    ast(obj.flag);
                    end_child(obj);
                    
            end
        end
    end
end