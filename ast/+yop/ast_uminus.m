classdef ast_uminus < yop.ast_unary_expression
    
    properties (Constant)
        name = 'uminus'
    end
    
    methods
        
        function obj = ast_uminus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
    end
    
end