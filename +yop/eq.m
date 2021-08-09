classdef eq < yop.binary_relation
    
    methods
        function obj = eq(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'eq(lhs, rhs)');
        end
    end
end