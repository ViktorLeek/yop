classdef ast_gt < yop.ast_binary_relation
    
    properties (Constant)
        name = 'gt'
    end
    
    methods
        function obj = ast_gt(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end