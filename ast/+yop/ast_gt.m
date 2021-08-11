classdef ast_gt < yop.ast_binary_relation
    
    methods
        function obj = ast_gt(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function ast(obj)
            xast(obj, 'gt(lhs, rhs)');
        end
    end
end