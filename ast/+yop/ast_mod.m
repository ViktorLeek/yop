classdef ast_mod < yop.ast_expression
    properties
        a
        m
    end
    methods
        function obj = ast_mod(a, m)
            obj@yop.ast_expression();
            obj.a = a;
            obj.m = m;
            obj.dim = size(mod(ones(size(a)), ones(size(m))));
        end
        
        function value = evaluate(obj)
            value = mod(evaluate(obj.a), evaluate(obj.m));
        end
        
        function draw(obj)
            fprintf('mod(a, m)\n');
            
            begin_child(obj);
            draw(obj.a);
            end_child(obj);
            
            last_child(obj);
            draw(obj.m);
            end_child(obj);
        end
    end
end