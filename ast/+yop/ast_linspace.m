classdef ast_linspace < yop.ast_expression
    properties
        x1
        x2
        n
    end
    methods
        function obj = ast_linspace(x1, x2, n)
            obj@yop.ast_expression();
            obj.x1 = x1;
            obj.x2 = x2;
            obj.n = n;
            obj.dim = [1, n];
        end
        
        function value = evaluate(obj)
            value = linspace(...
                evaluate(obj.x1), ...
                evaluate(obj.x2), ...
                evaluate(obj.n) ...
                );
        end
        
        function draw(obj)
            fprintf('linspace(x1, x2, n)\n');
            
            begin_child(obj);
            draw(obj.x1);
            end_child(obj);
            
            begin_child(obj);
            draw(obj.x2);
            end_child(obj);
            
            last_child(obj);
            draw(obj.n);
            end_child(obj);
        end
    end
end