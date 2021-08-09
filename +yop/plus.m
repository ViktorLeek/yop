classdef plus < yop.binary_expression
    
    methods
        function obj = plus(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'plus(lhs, rhs)');
        end
    end
end