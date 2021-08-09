classdef minus < yop.binary_expression
    
    methods
        function obj = minus(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'minus(lhs, rhs)');
        end
    end
end