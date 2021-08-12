classdef ast_le < yop.ast_binary_relation
    
    properties (Constant)
        name = 'le'
    end
    
    methods
        function obj = ast_le(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end