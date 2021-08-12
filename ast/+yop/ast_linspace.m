classdef ast_linspace < yop.ast_node
    properties
        x1
        x2
        n
    end
    methods
        function obj = ast_linspace(x1, x2, n)
            obj.x1 = x1;
            obj.x2 = x2;
            obj.n = n;
            obj.dim = [1, n];
        end
        function ast(obj)
            fprintf('linspace(x1, x2, n)\n');
            
            begin_child(obj);
            ast(obj.x1);
            end_child(obj);
            
            begin_child(obj);
            ast(obj.x2);
            end_child(obj);
            
            last_child(obj);
            ast(obj.n);
            end_child(obj);
        end
    end
end