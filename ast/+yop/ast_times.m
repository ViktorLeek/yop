classdef ast_times < yop.ast_binary_expression
    
    properties (Constant)
        name = 'times'
    end
    
    methods
        function obj = ast_times(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( times(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
end