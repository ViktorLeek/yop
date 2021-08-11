classdef ast_uplus < yop.ast_unary_expression
    
    methods
        
        function obj = ast_uplus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
    end
    
    methods % Printing
        
        function ast(obj)
            xast(obj, 'uplus(expr)');
        end
        
    end
    
end