classdef ast_alg < yop.ast_node
    properties
        var
    end
    methods
        function obj = ast_alg(var)
            obj.var = var;
            obj.dim = size(var);
        end
        function ast(obj)
            fprintf('alg(var)\n');
            last_child(obj);
            ast(obj.var);
            end_child(obj);
        end
    end
end