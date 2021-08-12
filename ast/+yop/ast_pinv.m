classdef ast_pinv < yop.ast_node
    properties
        A
    end
    methods
        function obj = ast_pinv(A)
            obj.A = A;
            obj.dim = size(pinv(ones(size(A))));
        end
        function ast(obj)
            fprintf('pinv(A)\n');
            last_child(obj);
            ast(obj.A);
            end_child(obj);
        end
    end
end