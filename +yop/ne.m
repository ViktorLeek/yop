classdef ne < yop.binary_relation
    
    methods
        function obj = ne(lhs, rhs)
            obj@yop.binary_relation(lhs, rhs);
        end
    end
    
    methods % Printing
        function print(obj)
            xprint(obj, 'ne(lhs, rhs)');
        end
    end
end