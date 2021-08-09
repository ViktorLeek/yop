classdef uminus < yop.unary_expression
    
    methods
        
        function obj = uminus(expr)
            obj@yop.unary_expression(expr);
        end
        
    end
    
    methods % Printing
        
        function print(obj)
            xprint(obj, 'uminus(expr)');
        end
        
    end
    
end