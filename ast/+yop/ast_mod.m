classdef ast_mod < yop.ast_node
    properties
        a
        m
    end
    methods
        function obj = ast_mod(a, m)
            obj.a = a;
            obj.m = m;
            obj.dim = size(mod(ones(size(a)), ones(size(m))));
        end
        function ast(obj)
            fprintf('mod(a, m)\n');
            
            begin_child(obj);
            ast(obj.a);
            end_child(obj);
            
            last_child(obj);
            ast(obj.m);
            end_child(obj);
        end
    end
end