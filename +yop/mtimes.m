classdef mtimes < yop.binary_expression
    
    methods
        function obj = mtimes(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
            obj.dim = size( mtimes(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'mtimes(lhs, rhs)');
        end
    end
end