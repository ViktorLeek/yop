classdef ast_pinv < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_pinv(A)
            obj.A = A;
            obj.dim = size(pinv(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = pinv(evaluate(obj.A));
        end
        
        function draw(obj)
            fprintf('pinv(A)\n');
            last_child(obj);
            draw(obj.A);
            end_child(obj);
        end
    end
end