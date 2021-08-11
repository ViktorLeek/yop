classdef ast_ne < yop.ast_binary_relation
    
    methods
        function obj = ast_ne(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'ne(lhs, rhs)');
        end
    end
end