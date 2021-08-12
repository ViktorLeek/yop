classdef ast_eq < yop.ast_binary_relation
    
    properties (Constant)
        name = 'eq'
    end
    
    methods
        function obj = ast_eq(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end