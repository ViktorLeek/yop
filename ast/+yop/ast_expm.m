classdef ast_expm < yop.ast_node
    properties
        A
    end
    methods
        function obj = ast_expm(A)
            obj.A = A;
            obj.dim = size(expm(ones(size(A))));
        end
        
        function value = evaluate(obj)
            value = expm(evaluate(obj.A));
        end
        
        function ast(obj)
            fprintf('expm(A)\n');
            last_child(obj);
            ast(obj.A);
            end_child(obj);
        end
    end
end