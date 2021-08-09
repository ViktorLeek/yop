classdef plus < yop.binary_expression
    
    methods
        function obj = plus(lhs, rhs)
            obj@yop.binary_expression(lhs, rhs);
            obj.dim = size( plus(ones(size(lhs)), ones(size(rhs))) );
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'plus(lhs, rhs)');
        end
    end
end