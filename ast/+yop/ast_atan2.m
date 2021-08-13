classdef ast_atan2 < yop.ast_node
    properties
        y
        x
    end
    methods
        function obj = ast_atan2(y, x)
            obj.y = y;
            obj.x = x;
            obj.dim = size(atan2(ones(size(y)), ones(size(x))));
        end
        
        function value = evaluate(obj)
            value = atan2(evaluate(obj.y), evaluate(obj.x));
        end
        
        function ast(obj)
            fprintf('atan2(y, x)\n');
            
            begin_child(obj);
            ast(obj.y);
            end_child(obj);
            
            last_child(obj);
            ast(obj.x);
            end_child(obj);
        end
    end
end