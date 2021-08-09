classdef mpower < yop.binary_expression
    
    methods
        function obj = mpower(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'mpower(lhs, rhs)');
        end
    end
end