classdef mldivide < yop.binary_expression
    
    methods
        function obj = mldivide(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'mldivide(lhs, rhs)');
        end
    end
end