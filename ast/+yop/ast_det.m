classdef ast_det < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_det(A)
            obj.A = A;
            obj.dim = size(det(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = det(evaluate(obj.A));
        end
        
        function draw(obj)
            fprintf('det(A)\n');
            last_child(obj);
            draw(obj.A);
            end_child(obj);
        end
    end
end