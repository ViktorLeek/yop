classdef mrdivide < yop.binary_expression
    
    methods
        function obj = mrdivide(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
            obj.dim = size( mrdivide(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'mrdivide(lhs, rhs)');
        end
    end
end