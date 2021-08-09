classdef ge < yop.binary_relation
    
    methods
        function obj = ge(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'ge(lhs, rhs)');
        end
    end
end