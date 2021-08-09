classdef uplus < yop.unary_expression
    
    methods
        
        function obj = uplus(expr)
            obj@yop.unary_expression(expr);
            obj.dim = size(expr);
        end
        
    end
    
    methods % Printing
        
        function print(obj)
            xprint(obj, 'uplus(expr)');
        end
        
    end
    
end