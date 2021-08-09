classdef gt < yop.binary_relation
    
    methods
        function obj = gt(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'gt(lhs, rhs)');
        end
    end
end