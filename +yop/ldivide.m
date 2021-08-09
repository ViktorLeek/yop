classdef ldivide < yop.binary_expression
    
    methods
        function obj = ldivide(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
            obj.dim = size( ldivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'ldivide(lhs, rhs)');
        end
    end
end