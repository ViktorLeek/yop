classdef power < yop.binary_expression
    
    methods
        function obj = power(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'power(lhs, rhs)');
        end
    end
end