classdef le < yop.binary_relation
    
    methods
        function obj = le(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'le(lhs, rhs)');
        end
    end
end