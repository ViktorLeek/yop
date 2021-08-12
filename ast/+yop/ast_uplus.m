classdef ast_uplus < yop.ast_unary_expression
    
    properties (Constant)
        name = 'uplus'
    end
    
    methods
        
        function obj = ast_uplus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
    end
    
end