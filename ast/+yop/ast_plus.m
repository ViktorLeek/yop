classdef ast_plus < yop.ast_binary_expression
    
    properties (Constant)
        name = 'plus'
    end
    
    methods
        function obj = ast_plus(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( plus(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
end