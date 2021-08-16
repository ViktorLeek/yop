classdef ast_cumsum < yop.ast_expression
    properties
        A
        d
        nargs
    end
    methods
        function obj = ast_dot(A, d)
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
    end
end