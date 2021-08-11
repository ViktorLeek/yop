classdef ast_ge < yop.ast_binary_relation
    
    methods
        function obj = ast_ge(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'ge(lhs, rhs)');
        end
    end
end