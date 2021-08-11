classdef ast_eq < yop.ast_binary_relation
    
    methods
        function obj = ast_eq(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'eq(lhs, rhs)');
        end
    end
end