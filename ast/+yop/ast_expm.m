classdef ast_expm < yop.ast_expression
    properties
        A
    end
    methods
        function obj = ast_expm(A)
            obj@yop.ast_expression();
            obj.A = A;
            obj.dim = size(expm(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = expm(evaluate(obj.A));
        end
        
        function draw(obj)
            fprintf('expm(A)\n');
            last_child(obj);
            draw(obj.A);
            end_child(obj);
        end
    end
end