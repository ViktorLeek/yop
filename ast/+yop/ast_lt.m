classdef ast_lt < yop.ast_binary_relation
    
    properties (Constant)
        name = 'lt'
    end
    
    methods
        function obj = ast_lt(lhs, rhs)
            obj@yop.ast_binary_relation(lhs, rhs);
        end
        
        function value = evaluate(obj)
            value = lt(evaluate(obj.lhs), evaluate(obj.rhs));
        end
    end
    
end