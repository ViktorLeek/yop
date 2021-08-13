classdef ast_det < yop.ast_node
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
        
        function ast(obj)
            fprintf('det(A)\n');
            last_child(obj);
            ast(obj.A);
            end_child(obj);
        end
    end
end