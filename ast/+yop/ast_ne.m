classdef ast_ne < yop.ast_binary_relation
    
    properties (Constant)
        name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end