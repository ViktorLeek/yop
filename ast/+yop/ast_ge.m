classdef ast_ge < yop.ast_binary_relation
    
    properties (Constant)
        name = 'ge'
    end
    
    methods
        function obj = ast_ge(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end