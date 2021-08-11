classdef ast_lt < yop.ast_binary_relation
    
    methods
        function obj = ast_lt(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'lt(lhs, rhs)');
        end
    end
end