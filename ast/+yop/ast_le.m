classdef ast_le < yop.ast_binary_relation
    
    methods
        function obj = ast_le(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'le(lhs, rhs)');
        end
    end
end