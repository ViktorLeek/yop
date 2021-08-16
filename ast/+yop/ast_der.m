classdef ast_der < yop.ast_expression
    properties
        var
    end
    methods
        function obj = ast_der(var)
            obj.var = var;
            obj.dim = size(var);
        end
        
        function value = evaluate(obj)
            value = evaluate(obj.var);
        end
        
        function draw(obj)
            fprintf('der(var)\n');
            last_child(obj);
            draw(obj.var);
            end_child(obj);
        end
    end
end