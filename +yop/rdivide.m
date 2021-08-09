classdef rdivide < yop.binary_expression
    
    methods
        function obj = rdivide(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'rdivide(lhs, rhs)');
        end
    end
end