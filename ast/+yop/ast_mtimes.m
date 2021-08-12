classdef ast_mtimes < yop.ast_binary_expression
    
    properties (Constant)
        name = 'mtimes'
    end
    
    methods
        function obj = ast_mtimes(lhs, rhs)
            obj@yop.ast_binary_expression(lhs, rhs);
            obj.dim = size( mtimes(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
end