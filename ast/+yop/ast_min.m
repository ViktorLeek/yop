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
            obj.A = A;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(min(ones(size(A))));
                    
                case 2
                    obj.B = B;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp));
                    
                case 3
                    obj.B = B;
                    obj.d = d;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(min(ones(size(A)), tmp, d));
                    
                case 4
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
        
        function ast(obj)
            switch obj.nargs
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