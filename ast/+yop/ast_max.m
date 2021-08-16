classdef ast_max < yop.ast_expression
    properties
        A
        B
        d
        flag
        nargs
    end
    methods
        function obj = ast_max(A, B, d, flag)
            obj.A = A;
            obj.nargs = nargin;
            switch nargin
                case 1
                    obj.dim = size(max(ones(size(A))));
                    
                case 2
                    obj.B = B;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(max(ones(size(A)), tmp));
                    
                case 3
                    obj.B = B;
                    obj.d = d;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(max(ones(size(A)), tmp, d));
                    
                case 4
                    obj.B = B;
                    obj.d = d;
                    obj.flag = flag;
                    if isa(B, 'yop.ast_node')
                        tmp = ones(size(B));
                    else
                        tmp = B;
                    end
                    obj.dim = size(max(ones(size(A)), tmp, d, flag));     
            end
            
        end
        
        function value = evaluate(obj)
            switch obj.nargs
                case 1
                    value = max(evaluate(obj.A));
                    
                case 2
                    value = max(evaluate(obj.A), evaluate(obj.B));
                    
                case 3
                    value = max(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d) ...
                        );
                    
                case 4
                    value = max(...
                        evaluate(obj.A), ...
                        evaluate(obj.B), ...
                        evaluate(obj.d), ...
                        evaluate(obj.flag) ...
                        );
            end
        end
        
        function draw(obj)
            switch obj.nargs
                case 1
                    fprintf('max(A)\n');
                    last_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                case 2
                    fprintf('max(A, B)\n');
                    
                    begin_child(obj);
                    draw(obj.A);
                    end_child(obj);
                    
                    last_child(obj);
                    draw(obj.B);
                    end_child(obj);
                    
                case 3
                    fprintf('max(A, B, dim)\n');
                    
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
                    fprintf('max(A, B, dim, flag)\n');
                    
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
    end
end