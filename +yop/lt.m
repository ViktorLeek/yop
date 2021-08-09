classdef lt < yop.binary_relation
    
    methods
        function obj = lt(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'lt(lhs, rhs)');
        end
    end
end