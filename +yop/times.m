classdef times < yop.binary_expression
    
    methods
        function obj = times(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'times(lhs, rhs)');
        end
    end
end