classdef ast_lt < yop.ast_binary_relation
    
    properties (Constant)
        name = 'lt'
    end
    
    methods
        function obj = ast_lt(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
    end
    
end