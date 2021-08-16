classdef ast_ne < yop.ast_relation
    
    properties (Constant)
        name = 'ne'
    end
    
    methods
        function obj = ast_ne(lhs, rhs)
            obj@yop.ast_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = ne(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end