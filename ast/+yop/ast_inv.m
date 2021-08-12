classdef ast_inv < yop.ast_node
    properties
        A
    end
    methods
        function obj = ast_inv(A)
            obj.A = A;
            obj.dim = size(inv(ones(size(A))));
        end
        function ast(obj)
            fprintf('inv(A)\n');
            last_child(obj);
            ast(obj.A);
            end_child(obj);
        end
    end
end