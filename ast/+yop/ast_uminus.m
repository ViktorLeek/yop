classdef ast_uminus < yop.ast_unary_expression
    
    methods
        
        function obj = ast_uminus(expr)
            obj@yop.ast_unary_expression(expr);
            obj.dim = size(expr);
        end
        
    end
    
    methods % Printing
        
        function ast(obj)
            xast(obj, 'uminus(expr)');
        end
        
    end
    
end