classdef ast_inv < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_inv(A)
            obj.A = A;
            obj.dim = size(inv(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = inv(evaluate(obj.A));
        end
        
        function draw(obj)
            fprintf('inv(A)\n');
            last_child(obj);
            draw(obj.A);
            end_child(obj);
        end
    end
end