classdef ast_der < yop.ast_node
    properties
        var
    end
    methods
        function obj = ast_der(var)
            obj.var = var;
            obj.dim = size(var);
        end
        function ast(obj)
            fprintf('der(var)\n');
            last_child(obj);
            ast(obj.var);
            end_child(obj);
        end
    end
end